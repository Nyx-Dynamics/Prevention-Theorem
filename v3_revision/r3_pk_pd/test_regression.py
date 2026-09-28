"""
R3 Regression Test (kinetically-aware refactor)
===============================================

Compares the v2 step-function PEP framework (eps_max=0.95, eps_mid=0.50,
eps_min=0 as constants) against the R3 kinetically-aware framework
(window-averaged eps_drug(t_PEP) * p_clearable(stage)) using the v2
multiscale model realizations.

Reports per route:
    - t_crit at eta=0.05 under v2 vs R3 across adherence
    - E_PEP at canonical timepoints
    - Envelope bound under Eq. 4
    - Hartford-specific evaluation (parenteral structural delay)

Author: A.C. Demidont, DO, AAHIVS (Nyx Dynamics LLC)
"""

from pathlib import Path
from math import log, erf, sqrt
import numpy as np
import pandas as pd

from effective_epsilon import (
    compute_eps_drug_curve, compute_E_PEP_r3_curve,
    T_CRIT_PARENTERAL_H, T_CRIT_MUCOSAL_H,
    DEFAULT_STAGE_CLEARABILITY,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
REALIZATIONS_CSV = (
    REPO_ROOT / 'SRC' / 'multiscale_model' / 'results_v3' /
    'heterogeneity_realizations.csv'
)
SUMMARY_CSV = (
    REPO_ROOT / 'SRC' / 'multiscale_model' / 'results_v3' /
    'heterogeneity_summary.csv'
)
ROUTE_FROM_V0 = {1: 'mucosal', 1000: 'parenteral'}


def load_realizations(V0: int, cv: float) -> pd.DataFrame:
    df = pd.read_csv(REALIZATIONS_CSV)
    sub = df[(df['V0'] == V0) & (df['cv'] == cv)]
    return sub.reset_index(drop=True)


def compute_seed_int_curves(
    df: pd.DataFrame, t_grid_h: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Empirical CDFs conditional on integration (v2 convention).
    Spontaneous extinction is handled separately via F_access in Eq. 4.
    """
    sub = df[df['extincted'] == False].copy()
    denom = len(sub)
    T_seed = sub['handoff_time_hours'].values
    T_int = sub['T_int_hours'].values
    P_seed = np.array([
        float(np.sum((T_seed <= t) & ~np.isnan(T_seed))) / denom
        for t in t_grid_h
    ])
    P_int = np.array([
        float(np.sum((T_int <= t) & ~np.isnan(T_int))) / denom
        for t in t_grid_h
    ])
    return P_seed, P_int


def find_t_crit(t_grid_h: np.ndarray, E_PEP: np.ndarray,
                eta: float = 0.05) -> float:
    below = np.where(E_PEP < eta)[0]
    if len(below) == 0:
        return float('inf')
    return float(t_grid_h[below[0]])


def _F_access(t_crit_h: float, F_median_h: float, GSD: float) -> float:
    """Fraction initiating PEP within t_crit under a log-normal delay."""
    if t_crit_h == float('inf'):
        return 1.0
    z = (log(t_crit_h) - log(F_median_h)) / log(GSD)
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def envelope_bound(t_crit_h: float, eps_max_at_zero: float,
                   F_median_h: float = 96.0, GSD: float = 2.0,
                   eta: float = 0.05) -> tuple[float, float]:
    """UPPER-BOUND CONVENTION of Equation 4. Retained for reproduction.

        bound = F_access(t_crit) * eps_max + (1 - F_access) * eta

    This reproduces the previously reported 11.3% at the canonical
    operating point and is kept unchanged so the manuscript-reproduction
    pathway is byte-identical.

    Note on the second term.  `eta` is the efficacy threshold that DEFINES
    t_crit (the level at which E_PEP is declared to have closed the
    window).  Here it is also serving as the residual efficacy credited to
    the (1 - F_access) fraction who never acquire drug inside the window.
    Under the v4 gating-event framing that fraction receives no drug, so
    its contribution is zero, and this expression is therefore the UPPER
    end of a range rather than a point estimate.  Use
    envelope_bound_range() for both ends; see numerical_claims_v4.csv
    entry `bound_envelope_par_range`.

    Under R3, eps_max here is eps_drug(t_PEP=0) — the PEP efficacy
    at t_PEP just after exposure, which is the upper-bound PEP-side
    efficacy term in Eq. 4. The biology p_clearable(0) = 1 means
    this equals eps_drug(0).
    """
    if t_crit_h == float('inf'):
        return 1.0, eps_max_at_zero
    F_access = _F_access(t_crit_h, F_median_h, GSD)
    bound = F_access * eps_max_at_zero + (1.0 - F_access) * eta
    return F_access, bound


def envelope_bound_range(t_crit_h: float, eps_max_at_zero: float,
                         F_median_h: float = 96.0, GSD: float = 2.0,
                         eta: float = 0.05) -> tuple[float, float, float]:
    """Distribution-free Equation 4 bound as a RANGE, not a point.

    Returns (F_access, lower, upper):

        lower = F_access(W) * eps_max
                  no residual protection after the window (v4 gating rule:
                  no drug acquired before t_crit implies E_PEP = 0)

        upper = F_access(W) * eps_max + (1 - F_access(W)) * eta
                  residual efficacy credited up to eta for the fraction
                  that misses the window; equals envelope_bound()

    The true bound lies between these, because someone initiating just
    past t_crit receives some efficacy — below eta by construction, but
    above zero.  At the canonical operating point (W = 34.5 h,
    eps_max = 0.95, F_median = 96 h, GSD = 2.0) this returns
    6.6% - 11.3%, which is the value reported in the CROI abstract.
    """
    if t_crit_h == float('inf'):
        return 1.0, eps_max_at_zero, eps_max_at_zero
    F_access = _F_access(t_crit_h, F_median_h, GSD)
    lower = F_access * eps_max_at_zero
    upper = lower + (1.0 - F_access) * eta
    return F_access, lower, upper


def envelope_band(windows_h, eps_max_at_zero: float = 0.95,
                  F_median_h: float = 96.0, GSD: float = 2.0,
                  eta: float = 0.05) -> pd.DataFrame:
    """Window-sensitivity band from envelope_bound_range().

    Generates the CROI graphic: one row per assumed PEP window W, with
    F_access(W) and the lower/upper ends of the distribution-free bound.
    """
    rows = []
    for W in windows_h:
        F, lo, hi = envelope_bound_range(W, eps_max_at_zero, F_median_h,
                                         GSD, eta)
        rows.append({'window_h': W, 'F_access': F,
                     'ceiling_lower_pct': lo * 100.0,
                     'ceiling_upper_pct': hi * 100.0})
    return pd.DataFrame(rows)


def v2_E_PEP(P_seed, P_int):
    """v2 step-function formula with hand-set constants."""
    return (1 - P_seed) * 0.95 + (P_seed - P_int) * 0.50 + P_int * 0.0


def run(V0: int, cv: float, t_crit_route_h: float):
    route = ROUTE_FROM_V0[V0]
    print(f"\n{'='*78}")
    print(f"R3 Kinetically-Aware Regression: {route} V0={V0} CV={cv} "
          f"(t_crit_route={t_crit_route_h}h)")
    print('='*78)

    df = load_realizations(V0, cv)
    summary = pd.read_csv(SUMMARY_CSV)
    can = summary[(summary['route'] == route) &
                  (summary['V0'] == V0) &
                  (summary['cv'] == cv)].iloc[0]
    v2_t_crit = can['t_crit_eta_05']
    p_extinct = can['p_extinct']
    n_int = can['n_int']

    print(f"\nMultiscale model state (commit 37e27ea):")
    print(f"  N_realizations: {len(df)}, P(extinct)={p_extinct}, "
          f"N_integrated={n_int}")
    print(f"  v2 canonical t_crit @ eta=0.05: {v2_t_crit:.1f} h")

    # t_PEP grid out to slightly past t_crit_route
    t_pep_grid = np.arange(0.0, t_crit_route_h + 5.0, 0.5)
    P_seed, P_int = compute_seed_int_curves(df, t_pep_grid)

    # v2 reference
    E_PEP_v2 = v2_E_PEP(P_seed, P_int)
    t_crit_v2 = find_t_crit(t_pep_grid, E_PEP_v2)

    # R3 kinetically-aware across adherence
    print(f"\n  {'Adh':>5} {'eps_drug(0)':>13} {'eps_drug(t_crit-1h)':>22}"
          f" {'t_crit_R3':>10} {'shift':>8} {'F_acc':>7} {'Bound':>7}")
    print('-'*78)

    # v2 row first
    F_v2, bnd_v2 = envelope_bound(t_crit_v2, 0.95)
    print(f"  {'v2':>5} {0.95:>13.4f} {'--':>22} "
          f"{t_crit_v2:>9.1f}h {0.0:>+7.1f}h {F_v2:>7.4f} {bnd_v2:>7.4f}")

    rows = []
    for adh in (1.00, 0.95, 0.90, 0.80, 0.70, 0.50, 0.30):
        r = compute_E_PEP_r3_curve(
            t_pep_grid, P_seed, P_int, t_crit_route_h,
            adherence=adh, n_replicates=20,
        )
        t_crit_r3 = find_t_crit(t_pep_grid, r['E_PEP'])
        eps_at_zero = float(r['eps_drug'][0])
        # eps at 1h before route window (for parenteral, t_PEP=33.5h)
        ref_idx = int(np.argmin(np.abs(t_pep_grid - (t_crit_route_h - 1.0))))
        eps_at_ref = float(r['eps_drug'][ref_idx])
        shift = t_crit_r3 - t_crit_v2
        F_r3, bnd_r3 = envelope_bound(t_crit_r3, eps_at_zero)
        rows.append({
            'adherence': adh, 'eps_drug_at_0': eps_at_zero,
            'eps_drug_late': eps_at_ref, 't_crit_R3_h': t_crit_r3,
            'shift_vs_v2_h': shift, 'F_access': F_r3, 'bound': bnd_r3,
        })
        print(f"  {adh:>5.2f} {eps_at_zero:>13.4f} {eps_at_ref:>22.4f} "
              f"{t_crit_r3:>9.1f}h {shift:>+7.1f}h "
              f"{F_r3:>7.4f} {bnd_r3:>7.4f}")

    # E_PEP at canonical timepoints under perfect adherence
    r_perf = compute_E_PEP_r3_curve(
        t_pep_grid, P_seed, P_int, t_crit_route_h,
        adherence=1.0, n_replicates=1,
    )
    print(f"\nE_PEP at canonical timepoints (perfect adherence vs v2):")
    print(f"  {'t (h)':>7} {'P_seed':>8} {'P_int':>8} "
          f"{'v2 E_PEP':>9} {'R3 E_PEP':>9} {'diff':>8}")
    for tq in (6.0, 12.0, 24.0, 34.5, 48.0, 60.5, 72.0):
        if tq > t_pep_grid[-1]:
            continue
        idx = int(np.argmin(np.abs(t_pep_grid - tq)))
        diff = r_perf['E_PEP'][idx] - E_PEP_v2[idx]
        print(f"  {tq:>7.1f} {P_seed[idx]:>8.4f} {P_int[idx]:>8.4f} "
              f"{E_PEP_v2[idx]:>9.4f} {r_perf['E_PEP'][idx]:>9.4f} "
              f"{diff:>+8.4f}")

    perf = rows[0]
    print(f"\n{'-'*78}")
    print(f"VERDICT  ({route})")
    print('-'*78)
    print(f"  v2 t_crit reproduced from realizations: {t_crit_v2:.1f}h "
          f"(canonical: {v2_t_crit:.1f}h)")
    print(f"  R3 t_crit at 100% adherence: {perf['t_crit_R3_h']:.1f}h "
          f"(shift {perf['shift_vs_v2_h']:+.1f}h)")
    print(f"  R3 t_crit at 70% adherence:  "
          f"{[r for r in rows if r['adherence']==0.70][0]['t_crit_R3_h']:.1f}h")
    print(f"  R3 t_crit at 50% adherence:  "
          f"{[r for r in rows if r['adherence']==0.50][0]['t_crit_R3_h']:.1f}h")
    print(f"  R3 envelope bound at perf adherence:  "
          f"{perf['bound']*100:.1f}%  (v2: {bnd_v2*100:.1f}%)")
    return rows


if __name__ == '__main__':
    par_rows = run(V0=1000, cv=0.3, t_crit_route_h=T_CRIT_PARENTERAL_H)
    muc_rows = run(V0=1, cv=0.3, t_crit_route_h=T_CRIT_MUCOSAL_H)

    # Hartford-specific check (parenteral, ~6h structural delay per
    # v2 supplement; the Hartford-specific delay is in the same band)
    print(f"\n\n{'='*78}")
    print("Hartford-specific evaluation (parenteral exposure)")
    print('='*78)
    print(f"Hartford expected efficacy under v2 was reported at 78.8% "
          f"at structural delay")
    print(f"Per v2 supplement: median-city structural delay ~ 6 h "
          f"(beta_1=1.2 h/%, beta_2=0.3 h/%)")
    print(f"Hartford-specific delay (LateDx 22.4%, Linkage 73.4%): "
          f"~6 + 1.2*(22.4-10) + 0.3*(90-73.4) = 6 + 14.88 + 4.98 = ~25.9 h")
    print(f"\nLet's compute R3 E_PEP at the Hartford structural delay "
          f"under perfect adherence:")
    t_pep_grid = np.arange(0.0, 40.0, 0.5)
    df_par = load_realizations(V0=1000, cv=0.3)
    P_seed, P_int = compute_seed_int_curves(df_par, t_pep_grid)
    r_perf_par = compute_E_PEP_r3_curve(
        t_pep_grid, P_seed, P_int, T_CRIT_PARENTERAL_H,
        adherence=1.0, n_replicates=1,
    )
    print(f"  {'delay (h)':>11} {'v2 E_PEP':>9} {'R3 E_PEP':>9}")
    for delay in (6.0, 10.0, 15.0, 20.0, 24.0, 25.9, 30.0):
        idx = int(np.argmin(np.abs(t_pep_grid - delay)))
        e_v2 = v2_E_PEP(P_seed, P_int)[idx]
        e_r3 = r_perf_par['E_PEP'][idx]
        marker = ' <-- Hartford' if abs(delay - 25.9) < 0.1 else ''
        print(f"  {delay:>11.1f} {e_v2:>9.4f} {e_r3:>9.4f}{marker}")
