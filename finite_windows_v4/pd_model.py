"""
PD Module for TDF/FTC/DTG PEP Regimen (window-aware refactor)
=============================================================

Maps drug concentrations to pharmacologic suppression of new HIV-1
infection events.

Under the kinetically-aware framing (post AC's clarification on
2026-05-26), the PD operates on the *normalized active fraction*
of steady-state — not on raw plasma concentration:

    For NRTIs (TDF -> TFV-DP, FTC -> FTC-TP):
        active(t) = intracellular fraction f(t)  (f_ss = 1.0)
    For INSTI (DTG):
        active(t) = C_plasma(t) / C_plasma_ss    (also normalized)

This puts all drugs on the same dimensionless scale, with EC50 also
expressed as a fraction of steady-state. The Hill function gives:

    eps_i(active) = active^n / (EC50_frac^n + active^n)

Combination across drugs via multiplicative survival:

    eps_combined = 1 - product_i (1 - eps_i)

EC50_frac calibration anchors:
    TFV-DP intracellular EC50_frac ~ 0.10
        (cellular IC50 ~100 fmol/10^6 PBMCs; SS ~150-200 fmol/10^6;
         Anderson 2011, Castillo-Mancilla 2016)
    FTC-TP intracellular EC50_frac ~ 0.10
        (similar order; Wang 2004, Anderson 2011)
    DTG plasma EC50_frac    ~ 0.036
        (= EC50 ~80 ng/mL / C_avg_ss ~2200 ng/mL;
         Min 2010, FDA Tivicay label)

Hill coefficient n = 1 (no cooperativity). Sensitivity to n in [1, 3]
available as an optional parameter.

Author: A.C. Demidont, DO, AAHIVS (Nyx Dynamics LLC)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np

from pk_model import DrugPK, TDF, FTC, DTG, concentration_timecourse


# ----------------------------------------------------------------------
# Per-drug PD parameter records
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class DrugPD:
    """Hill-function PD parameters on normalized active fraction.

    EC50_frac is the active-fraction value at which suppression = 0.5.
    At active = 1.0 (steady state), eps = 1/(EC50_frac + 1).
    """
    name: str
    EC50_frac: float
    hill_n: float = 1.0


# Default PD parameters paired with the PK drug records
PD_DEFAULTS = {
    'TDF (modeled as TFV active metabolite)': DrugPD(
        name='TFV-DP', EC50_frac=0.10, hill_n=1.0
    ),
    'FTC': DrugPD(
        name='FTC-TP', EC50_frac=0.10, hill_n=1.0
    ),
    'DTG': DrugPD(
        name='DTG', EC50_frac=0.036, hill_n=1.0
    ),
}


# ----------------------------------------------------------------------
# Hill function: single-drug suppression on normalized active
# ----------------------------------------------------------------------

def single_drug_suppression(
    active_fraction: np.ndarray,
    EC50_frac: float,
    hill_n: float = 1.0,
) -> np.ndarray:
    """Hill-function fraction of viral replication suppressed.

    eps(a) = a^n / (EC50^n + a^n)

    At a = EC50_frac, eps = 0.5. At a >> EC50_frac, eps -> 1.
    At a = 0, eps = 0.
    """
    a = np.maximum(active_fraction, 0.0)
    an = a ** hill_n
    EC50n = EC50_frac ** hill_n
    return an / (EC50n + an)


# ----------------------------------------------------------------------
# Multi-drug combination: multiplicative survival
# ----------------------------------------------------------------------

def combined_drug_suppression(
    eps_per_drug: Sequence[np.ndarray],
) -> np.ndarray:
    """Combine per-drug suppression via multiplicative survival.

    eps_combined = 1 - product_i (1 - eps_i)

    Independent molecular targets: NRTIs block RT; INSTI blocks
    integration. The combination model assumes drugs act on distinct
    steps so survival probabilities multiply.
    """
    survival = np.ones_like(eps_per_drug[0])
    for eps_i in eps_per_drug:
        survival = survival * (1.0 - eps_i)
    return 1.0 - survival


# ----------------------------------------------------------------------
# End-to-end: PK + PD trajectory
# ----------------------------------------------------------------------

def regimen_eps_timecourse(
    drugs_pk: Sequence[DrugPK] = (TDF, FTC, DTG),
    duration_h: float = 28 * 24,
    dose_interval_h: float = 24.0,
    sampling_resolution_h: float = 0.5,
    adherence: float = 1.0,
    rng_seed: int = 0,
    pd_overrides: dict | None = None,
) -> dict:
    """Compute the regimen's time-resolved suppression eps_drug(t).

    Returns a dict with:
        t_grid_h: time grid (hours from PEP initiation)
        eps_drug_t: combined regimen suppression at each time
        eps_per_drug: per-drug suppression trajectories
        active_per_drug: per-drug active-fraction trajectories
        drug_names: drug name list
    """
    pd_params = dict(PD_DEFAULTS)
    if pd_overrides:
        pd_params.update(pd_overrides)

    t_grid = None
    eps_per_drug = []
    active_per_drug = []

    for drug in drugs_pk:
        r = concentration_timecourse(
            drug,
            duration_h=duration_h,
            dose_interval_h=dose_interval_h,
            sampling_resolution_h=sampling_resolution_h,
            adherence=adherence,
            rng_seed=rng_seed,
        )
        if t_grid is None:
            t_grid = r['t_grid_h']

        pd_p = pd_params[drug.name]
        eps_i = single_drug_suppression(
            r['active'], pd_p.EC50_frac, pd_p.hill_n
        )

        eps_per_drug.append(eps_i)
        active_per_drug.append(r['active'])

    eps_drug_t = combined_drug_suppression(eps_per_drug)

    return {
        't_grid_h': t_grid,
        'eps_drug_t': eps_drug_t,
        'eps_per_drug': eps_per_drug,
        'active_per_drug': active_per_drug,
        'drug_names': [d.name for d in drugs_pk],
    }


# ----------------------------------------------------------------------
# Window-averaged eps_drug for given PEP-initiation delay
# ----------------------------------------------------------------------

def window_averaged_eps_drug(
    t_pep_h: float,
    t_crit_route_h: float,
    eps_timecourse: dict,
) -> float:
    """Average pharmacologic suppression delivered during the remaining
    route-specific integration window.

    For PEP initiated at delay t_PEP from exposure, the drug has the
    interval [t_PEP, t_crit_route] to act on integration prevention.
    This function returns the time-average of eps_combined over that
    interval, where time 0 in the eps_timecourse corresponds to the
    moment of PEP initiation:

        eps_drug(t_PEP) = (1 / Delta_t) * integral_0^Delta_t eps(tau) dtau
                         where Delta_t = t_crit_route - t_PEP

    If t_pep_h >= t_crit_route_h, no time remains and the function
    returns 0.
    """
    delta_t = t_crit_route_h - t_pep_h
    if delta_t <= 0:
        return 0.0

    t_grid = eps_timecourse['t_grid_h']
    eps = eps_timecourse['eps_drug_t']

    mask = t_grid <= delta_t
    if mask.sum() < 2:
        return 0.0

    t_sub = t_grid[mask]
    eps_sub = eps[mask]
    auc = float(np.trapz(eps_sub, t_sub))
    return auc / (t_sub[-1] - t_sub[0])


if __name__ == '__main__':
    # 1) Time-resolved eps_drug at perfect adherence
    print("PD Module Sanity Check (kinetically-aware framework)")
    print("=" * 70)
    result = regimen_eps_timecourse(adherence=1.0, rng_seed=0)
    t = result['t_grid_h']
    print(f"\neps_combined at canonical timepoints (perfect adherence):")
    print(f"  {'t (h)':>7} {'TFV':>7} {'FTC':>7} {'DTG':>7} {'combined':>9}")
    for tq in (2, 6, 12, 24, 34.5, 48, 60.5, 72, 96, 168):
        idx = int(np.argmin(np.abs(t - tq)))
        eps_tfv = result['eps_per_drug'][0][idx]
        eps_ftc = result['eps_per_drug'][1][idx]
        eps_dtg = result['eps_per_drug'][2][idx]
        eps_comb = result['eps_drug_t'][idx]
        print(f"  {tq:>7.1f} {eps_tfv:>7.3f} {eps_ftc:>7.3f} "
              f"{eps_dtg:>7.3f} {eps_comb:>9.4f}")

    # 2) Window-averaged eps_drug across PEP initiation delays
    print(f"\n\nWindow-averaged eps_drug(t_PEP) "
          f"-- key quantity for prevention efficacy")
    print("=" * 70)
    print(f"  {'t_PEP (h)':>10} {'parenteral (34.5h)':>20} "
          f"{'mucosal (60.5h)':>18}")
    print("  " + "-" * 56)
    for t_pep in (0, 6, 12, 18, 24, 30, 34, 40, 48, 55, 60):
        eps_par = window_averaged_eps_drug(t_pep, 34.5, result)
        eps_muc = window_averaged_eps_drug(t_pep, 60.5, result)
        par_str = f"{eps_par:>20.4f}" if t_pep < 34.5 else f"{'(past window)':>20}"
        muc_str = f"{eps_muc:>18.4f}" if t_pep < 60.5 else f"{'(past window)':>18}"
        print(f"  {t_pep:>10.1f} {par_str} {muc_str}")
