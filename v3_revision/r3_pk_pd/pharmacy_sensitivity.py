"""
pharmacy_sensitivity.py — v4 envelope-corridor framing for pharmacy access

This replaces the v3 scalar "envelope bound" calculation
   bound = F_access(t_crit) * eps_max + (1 - F_access(t_crit)) * eps_min
which mixed a population-level access distribution with per-patient PK
efficacy and used a phantom eps_min = 0.05 floor inconsistent with the
gating-event causal chain.

Under v4, pharmacy access is modeled as an UPSTREAM GATING EVENT:
patients who do not acquire medication before t_crit never receive drug,
PK is irrelevant, and E_PEP = 0 by construction.  The "envelope" is a 2D
corridor in (t_acq, E_PEP) space, bounded above by perfect-adherence PK
and below by low-adherence PK.  The corridor is a fixed property of viral
kinetics + drug PK; pharmacy delays slide cities rightward along it.

Reuses the canonical multiscale loading and ECDF construction from
test_regression.py byte-identical (same realizations CSV, same P_seed and
P_int curves) so numerical results match the canonical PK pipeline.

Outputs (in v3_revision/results/pharmacy_sensitivity_corrected/):
    envelope_corridor.csv             — t_acq grid x (upper, lower) curves
    city_envelope_positions.csv       — long format, one row per (city, dt_pharm)
    pharmacy_sensitivity_results.csv  — legacy schema, perfect-adherence only
    pharmacy_displacement_summary.csv — N cities past t_crit at each dt_pharm

A note on Hartford's structural delay
-------------------------------------
This script reads each city's structural_delay_h directly from
city_pep_efficacy_results.csv (Hartford = 24.4h, matching the manuscript
prose).  The prior pharmacy_sensitivity.py used an internal reconstruction
that produced Hartford = 27.2h, inconsistent with the city CSV and the
manuscript.  Under the corrected (24.4h) value, Hartford crosses t_crit_R3
at dt_pharm ~10h rather than the prior reported ~8h.  This brings the
displacement story into alignment with the canonical structural delay
table used elsewhere in the paper.

Author: A. C. Demidont, DO, AAHIVS (Nyx Dynamics LLC)
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from effective_epsilon import (
    compute_E_PEP_r3_curve,
    T_CRIT_PARENTERAL_H,
    DEFAULT_STAGE_CLEARABILITY,  # noqa: F401  (imported for parity with test_regression.py)
)


# ---------------------------------------------------------------------------
# Paths (mirror test_regression.py)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
REALIZATIONS_CSV = (
    REPO_ROOT / 'SRC' / 'multiscale_model' / 'results_v3' /
    'heterogeneity_realizations.csv'
)
SUMMARY_CSV = (
    REPO_ROOT / 'SRC' / 'multiscale_model' / 'results_v3' /
    'heterogeneity_summary.csv'
)

# Candidate locations for the city panel CSV.  First match wins.
CITY_CSV_CANDIDATES = [
    REPO_ROOT / 'results' / 'city_analysis' / 'city_pep_efficacy_results.csv',
    REPO_ROOT / 'city_pep_efficacy_results.csv',
    REPO_ROOT / 'v3_revision' / 'results' / 'city_pep_efficacy_results.csv',
]

OUT_DIR = REPO_ROOT / 'v3_revision' / 'results' / 'pharmacy_sensitivity_corrected'
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Canonical constants (parenteral baseline)
# ---------------------------------------------------------------------------
ROUTE = 'parenteral'
V0_PARENTERAL = 1000
CV_PARENTERAL = 0.3
T_CRIT_BIOLOGICAL = T_CRIT_PARENTERAL_H  # 34.5 h

# Adherence levels defining the corridor.  Upper = perfect; lower = bottom
# of the realistic range AC swept in test_regression.py.
ADHERENCE_UPPER = 1.0
ADHERENCE_LOWER = 0.30

# t_acq grid: 0 to slightly past biological t_crit (matches test_regression.py).
T_ACQ_GRID = np.arange(0.0, T_CRIT_BIOLOGICAL + 5.0, 0.5)

# Pharmacy delay sweep (preserved from v3 for continuity).
PHARMACY_DELAYS_H = [0, 2, 4, 6, 8, 10, 12]

# E_PEP threshold defining the "off corridor" line.  Matches test_regression
# convention (find_t_crit at eta=0.05).
ETA = 0.05


# ---------------------------------------------------------------------------
# Multiscale loaders — verbatim from test_regression.py
# ---------------------------------------------------------------------------
def load_realizations(V0: int, cv: float) -> pd.DataFrame:
    df = pd.read_csv(REALIZATIONS_CSV)
    sub = df[(df['V0'] == V0) & (df['cv'] == cv)]
    return sub.reset_index(drop=True)


def compute_seed_int_curves(
    df: pd.DataFrame, t_grid_h: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Empirical CDFs conditional on integration (v2 convention).
    Spontaneous extinction is handled separately upstream.
    """
    sub = df[df['extincted'] == False].copy()  # noqa: E712  (pandas idiom)
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
                eta: float = ETA) -> float:
    below = np.where(E_PEP < eta)[0]
    if len(below) == 0:
        return float('inf')
    return float(t_grid_h[below[0]])


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_city_csv() -> Path:
    for p in CITY_CSV_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(
        "Could not locate city_pep_efficacy_results.csv.  Searched: "
        + ", ".join(str(p) for p in CITY_CSV_CANDIDATES)
    )


def interp_at(t_query: float, t_grid: np.ndarray, curve: np.ndarray) -> float:
    """Linear interpolation of a curve at an arbitrary t."""
    if t_query <= t_grid[0]:
        return float(curve[0])
    if t_query >= t_grid[-1]:
        return float(curve[-1])
    return float(np.interp(t_query, t_grid, curve))


# ---------------------------------------------------------------------------
# Core: build corridor + city positions
# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("S14 PHARMACY ACCESS — v4 envelope-corridor framing")
    print(f"Route: {ROUTE} (V0={V0_PARENTERAL}, CV={CV_PARENTERAL})")
    print("=" * 80)

    # ---- Multiscale realizations (same data path as test_regression.py)
    df_par = load_realizations(V0_PARENTERAL, CV_PARENTERAL)
    print(f"\nLoaded {len(df_par)} realizations from {REALIZATIONS_CSV.name}")

    P_seed, P_int = compute_seed_int_curves(df_par, T_ACQ_GRID)

    # ---- Compute upper and lower corridor curves
    print(f"\nComputing envelope corridor:")
    print(f"  upper curve at rho={ADHERENCE_UPPER}")
    upper = compute_E_PEP_r3_curve(
        T_ACQ_GRID, P_seed, P_int, T_CRIT_BIOLOGICAL,
        adherence=ADHERENCE_UPPER, n_replicates=20,
    )
    print(f"  lower curve at rho={ADHERENCE_LOWER}")
    lower = compute_E_PEP_r3_curve(
        T_ACQ_GRID, P_seed, P_int, T_CRIT_BIOLOGICAL,
        adherence=ADHERENCE_LOWER, n_replicates=20,
    )

    t_crit_upper = find_t_crit(T_ACQ_GRID, upper['E_PEP'])
    t_crit_lower = find_t_crit(T_ACQ_GRID, lower['E_PEP'])
    print(f"\nAdherence-dependent t_crit (where E_PEP first drops below eta={ETA}):")
    print(f"  t_crit at rho={ADHERENCE_UPPER}: {t_crit_upper:.1f} h")
    print(f"  t_crit at rho={ADHERENCE_LOWER}: {t_crit_lower:.1f} h")

    # ---- Save corridor
    corridor_df = pd.DataFrame({
        't_acq_h': T_ACQ_GRID,
        'E_PEP_upper': upper['E_PEP'],
        'E_PEP_lower': lower['E_PEP'],
        'envelope_width': upper['E_PEP'] - lower['E_PEP'],
        'eps_drug_upper': upper['eps_drug'],
        'eps_drug_lower': lower['eps_drug'],
        'biology_term': upper['biology_term'],  # same for both adherence values
        'P_seed': P_seed,
        'P_int': P_int,
        'past_tcrit_upper': T_ACQ_GRID >= t_crit_upper,
        'past_tcrit_lower': T_ACQ_GRID >= t_crit_lower,
    })
    corridor_df.to_csv(OUT_DIR / 'envelope_corridor.csv', index=False)
    print(f"\nSaved {OUT_DIR / 'envelope_corridor.csv'}  ({len(corridor_df)} rows)")

    # ---- Load city panel
    city_csv = find_city_csv()
    cities_df = pd.read_csv(city_csv)
    print(f"\nLoaded city panel from {city_csv} ({len(cities_df)} cities)")

    if 'structural_delay_h' not in cities_df.columns:
        raise ValueError(
            f"city CSV at {city_csv} must contain a 'structural_delay_h' column. "
            f"Got columns: {list(cities_df.columns)}"
        )

    # ---- Compute city positions on the corridor
    city_rows = []
    legacy_rows = []
    for _, city in cities_df.iterrows():
        dt_struct = float(city['structural_delay_h'])
        for dt_pharm in PHARMACY_DELAYS_H:
            t_acq = dt_struct + dt_pharm
            e_up = interp_at(t_acq, T_ACQ_GRID, upper['E_PEP'])
            e_lo = interp_at(t_acq, T_ACQ_GRID, lower['E_PEP'])
            past_tcrit = t_acq >= t_crit_upper

            city_rows.append({
                'city': city['city'],
                'state': city.get('state', ''),
                'delta_t_struct_h': dt_struct,
                'delta_t_pharm_h': dt_pharm,
                't_acq_h': t_acq,
                'E_PEP_upper': e_up,
                'E_PEP_lower': e_lo,
                'past_tcrit': past_tcrit,
            })

            # Legacy schema (perfect adherence only) for back-compat
            retention = (
                float(city['linkage_to_care_pct']) / 100.0
                if 'linkage_to_care_pct' in cities_df.columns
                else float('nan')
            )
            legacy_rows.append({
                'city': city['city'],
                'state': city.get('state', ''),
                'late_dx_pct': city.get('late_dx_pct', float('nan')),
                'retention': retention,
                'delta_t_struct_h': dt_struct,
                'delta_t_pharmacy_h': dt_pharm,
                'delta_t_total_h': t_acq,
                'E_PEP_R3': e_up,
                'past_t_crit': past_tcrit,
            })

    city_positions_df = pd.DataFrame(city_rows)
    city_positions_df.to_csv(OUT_DIR / 'city_envelope_positions.csv', index=False)
    print(f"Saved {OUT_DIR / 'city_envelope_positions.csv'}  ({len(city_positions_df)} rows)")

    legacy_df = pd.DataFrame(legacy_rows)
    legacy_df.to_csv(OUT_DIR / 'pharmacy_sensitivity_results.csv', index=False)
    print(f"Saved {OUT_DIR / 'pharmacy_sensitivity_results.csv'}  ({len(legacy_df)} rows)")

    # ---- Displacement summary (replaces the misleading scalar bound)
    summary_rows = []
    for dt_pharm in PHARMACY_DELAYS_H:
        subset = city_positions_df[city_positions_df['delta_t_pharm_h'] == dt_pharm]
        n_past = int(subset['past_tcrit'].sum())
        n_total = len(subset)
        displaced_names = sorted(subset.loc[subset['past_tcrit'], 'city'].tolist())
        summary_rows.append({
            'delta_t_pharm_h': dt_pharm,
            'n_cities_displaced_past_tcrit': n_past,
            'n_cities_total': n_total,
            'fraction_displaced': n_past / n_total if n_total else 0.0,
            'median_t_acq_h': float(subset['t_acq_h'].median()),
            'mean_E_PEP_upper': float(subset['E_PEP_upper'].mean()),
            'mean_E_PEP_lower': float(subset['E_PEP_lower'].mean()),
            'displaced_city_names': ', '.join(displaced_names) or '-',
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / 'pharmacy_displacement_summary.csv', index=False)
    print(f"Saved {OUT_DIR / 'pharmacy_displacement_summary.csv'}  ({len(summary_df)} rows)")

    # ---- Console summary
    print("\n" + "=" * 80)
    print("PHARMACY DISPLACEMENT SUMMARY  (corridor framing, supersedes v3 scalar)")
    print("=" * 80)
    print(summary_df.to_string(index=False))

    print("\nKey points:")
    print(f"  - The envelope corridor is bounded by E_PEP at rho=1.0 (upper) and")
    print(f"    rho=0.30 (lower), with adherence-dependent t_crit cliffs at")
    print(f"    {t_crit_upper:.1f}h and {t_crit_lower:.1f}h.")
    print(f"  - The corridor is a property of viral kinetics + drug PK and does")
    print(f"    NOT change when pharmacy delays are added.  Cities slide rightward.")
    print(f"  - Cities whose t_acq exceeds t_crit are off the corridor (E_PEP = 0).")

    # ---- Hartford-specific evaluation (parallel to test_regression.py)
    hartford = city_positions_df[city_positions_df['city'].str.contains(
        'Hartford', case=False, na=False
    )]
    if not hartford.empty:
        print("\nHartford trajectory (city CSV structural_delay_h):")
        print(f"  Delay base: {hartford['delta_t_struct_h'].iloc[0]:.1f} h")
        for _, row in hartford.iterrows():
            marker = '  <-- off corridor' if row['past_tcrit'] else ''
            print(f"    dt_pharm = {int(row['delta_t_pharm_h']):2d}h  "
                  f"t_acq = {row['t_acq_h']:5.1f}h  "
                  f"E_PEP_upper = {row['E_PEP_upper']:.4f}  "
                  f"E_PEP_lower = {row['E_PEP_lower']:.4f}{marker}")
        print(f"  (Note: t_crit_upper = {t_crit_upper:.1f}h; Hartford crosses at the "
              f"first dt_pharm where t_acq >= t_crit_upper.)")


if __name__ == '__main__':
    main()
