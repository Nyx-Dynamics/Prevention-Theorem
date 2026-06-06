"""
Effective Epsilon Module for R3 (kinetically-aware window-averaged)
===================================================================

Implements the kinetically-correct separation of pharmacology and
stage biology, per AC's framing (2026-05-26):

  - The biological integration window is route-specific and given by
    the multiscale model (t_crit_parenteral ~ 34.5 h, t_crit_mucosal
    ~ 60.5 h at eta=0.05, V0/CV canonical).
  - The pharmacologic action delivered by PEP depends on how rapidly
    drugs reach therapeutic intracellular concentration during the
    REMAINING window [t_PEP, t_crit_route] after PEP initiation.
  - Drug onset kinetics are deterministic at the population scale
    (no clinically meaningful inter-individual PK variability to
    model); use a single canonical onset trajectory per drug.

Combined formulation:

    eps(stage, t_PEP) = eps_drug(t_PEP) * p_clearable(stage)

    eps_drug(t_PEP) = (1 / Delta_t) * integral_0^Delta_t eps_combined(tau) dtau
                      where Delta_t = t_crit_route - t_PEP

This puts the kinetic onset competition at the center of the model:
late PEP initiation leaves insufficient time for drug accumulation
inside the integration window, dropping eps_drug toward zero as
t_PEP approaches t_crit_route.

Calibration of p_clearable:
    p_clearable(Z=0) = 1.000   (no commitments at pre-seeding)
    p_clearable(Z=1) = 0.502   (recovers v2 eps_mid = 0.50 at the
                                hypothetical perfect-PK ceiling)
    p_clearable(Z=2) = 0.000   (post-integration irreversibility)

Author: A.C. Demidont, DO, AAHIVS (Nyx Dynamics LLC)
"""

from dataclasses import dataclass
from typing import Sequence, Optional
import numpy as np

from pk_model import DrugPK, TDF, FTC, DTG
from pd_model import (
    regimen_eps_timecourse, window_averaged_eps_drug, DrugPD, PD_DEFAULTS,
)


# Route-specific biological integration window from v2 canonical
# numerical_claims_v2.csv (multiscale model, commit 37e27ea):
T_CRIT_PARENTERAL_H = 34.5      # V0=10^3, CV=0.3, eta=0.05
T_CRIT_MUCOSAL_H = 60.5         # V0=1,    CV=0.3, eta=0.05


@dataclass(frozen=True)
class StageClearability:
    """Probability that the trajectory is still rescuable at stage Z."""
    p_Z0_pre_seeding: float = 1.000
    p_Z1_seeded_pre_int: float = 0.502
    p_Z2_integrated: float = 0.000

    def p(self, stage: int) -> float:
        if stage == 0:
            return self.p_Z0_pre_seeding
        elif stage == 1:
            return self.p_Z1_seeded_pre_int
        elif stage == 2:
            return self.p_Z2_integrated
        raise ValueError(f"stage must be 0, 1, or 2; got {stage}")


DEFAULT_STAGE_CLEARABILITY = StageClearability()


# ----------------------------------------------------------------------
# Canonical pharmacologic trajectory and window-averaged eps_drug
# ----------------------------------------------------------------------

def compute_eps_drug_curve(
    t_pep_grid_h: np.ndarray,
    t_crit_route_h: float,
    adherence: float = 1.0,
    drugs_pk: Sequence[DrugPK] = (TDF, FTC, DTG),
    rng_seed: int = 0,
    pd_overrides: Optional[dict] = None,
    n_replicates: int = 20,
) -> np.ndarray:
    """eps_drug(t_PEP) for each t_PEP on the grid.

    At adherence=1.0 the result is deterministic (single replicate).
    For adherence<1.0 we average over n_replicates missed-dose draws
    to reduce stochastic noise.
    """
    if adherence >= 1.0 - 1e-9:
        n_replicates = 1

    eps_curves = []
    for s in range(n_replicates):
        timecourse = regimen_eps_timecourse(
            drugs_pk=drugs_pk,
            adherence=adherence,
            rng_seed=rng_seed + s,
            pd_overrides=pd_overrides,
        )
        curve = np.array([
            window_averaged_eps_drug(t_pep, t_crit_route_h, timecourse)
            for t_pep in t_pep_grid_h
        ])
        eps_curves.append(curve)

    return np.mean(eps_curves, axis=0)


# ----------------------------------------------------------------------
# E_PEP curve under the R3 kinetic framework
# ----------------------------------------------------------------------

def compute_E_PEP_r3_curve(
    t_pep_grid_h: np.ndarray,
    P_seed: np.ndarray,
    P_int: np.ndarray,
    t_crit_route_h: float,
    adherence: float = 1.0,
    drugs_pk: Sequence[DrugPK] = (TDF, FTC, DTG),
    stage_clearability: StageClearability = DEFAULT_STAGE_CLEARABILITY,
    rng_seed: int = 0,
    pd_overrides: Optional[dict] = None,
    n_replicates: int = 20,
) -> dict:
    """E_PEP(t_PEP) under R3 kinetically-aware framework.

    E_PEP(t_PEP) = eps_drug(t_PEP) * [
                       (1 - P_seed) * p_clearable(0)
                     + (P_seed - P_int) * p_clearable(1)
                   ]

    P_int(t_PEP) * 0 from p_clearable(2) drops out.

    Returns dict with eps_drug curve and E_PEP curve.
    """
    eps_drug_curve = compute_eps_drug_curve(
        t_pep_grid_h, t_crit_route_h, adherence=adherence,
        drugs_pk=drugs_pk, rng_seed=rng_seed,
        pd_overrides=pd_overrides, n_replicates=n_replicates,
    )

    p0 = stage_clearability.p(0)
    p1 = stage_clearability.p(1)
    biology_term = (1.0 - P_seed) * p0 + (P_seed - P_int) * p1

    E_PEP = eps_drug_curve * biology_term

    return {
        't_pep_grid_h': t_pep_grid_h,
        'eps_drug': eps_drug_curve,
        'biology_term': biology_term,
        'E_PEP': E_PEP,
        'adherence': adherence,
        't_crit_route_h': t_crit_route_h,
    }


if __name__ == '__main__':
    print("Effective Epsilon (kinetically-aware) Sanity Check")
    print("=" * 70)

    t_pep_grid = np.arange(0.0, 80.0, 0.5)

    print(f"\neps_drug(t_PEP) for parenteral window (t_crit = "
          f"{T_CRIT_PARENTERAL_H} h):")
    print(f"  {'t_PEP':>8} {'100%':>9} {'90%':>9} {'70%':>9} {'50%':>9}")
    eps_drugs_par = {}
    for adh in (1.00, 0.90, 0.70, 0.50):
        eps_drugs_par[adh] = compute_eps_drug_curve(
            t_pep_grid, T_CRIT_PARENTERAL_H, adherence=adh, n_replicates=20
        )
    for t in (0, 6, 12, 18, 24, 30, 33):
        idx = int(np.argmin(np.abs(t_pep_grid - t)))
        row = [f"  {t:>8.1f}"]
        for adh in (1.00, 0.90, 0.70, 0.50):
            row.append(f" {eps_drugs_par[adh][idx]:>9.4f}")
        print(''.join(row))

    print(f"\neps_drug(t_PEP) for mucosal window (t_crit = "
          f"{T_CRIT_MUCOSAL_H} h):")
    print(f"  {'t_PEP':>8} {'100%':>9} {'90%':>9} {'70%':>9} {'50%':>9}")
    eps_drugs_muc = {}
    for adh in (1.00, 0.90, 0.70, 0.50):
        eps_drugs_muc[adh] = compute_eps_drug_curve(
            t_pep_grid, T_CRIT_MUCOSAL_H, adherence=adh, n_replicates=20
        )
    for t in (0, 12, 24, 36, 48, 55, 59):
        idx = int(np.argmin(np.abs(t_pep_grid - t)))
        row = [f"  {t:>8.1f}"]
        for adh in (1.00, 0.90, 0.70, 0.50):
            row.append(f" {eps_drugs_muc[adh][idx]:>9.4f}")
        print(''.join(row))
