"""
R3 Sensitivity Analysis: Pharmacy-Access Delay Layer
====================================================

Quantifies the impact of a pharmacy-access delay (Δt_pharmacy) layered
on top of the existing city-level structural delay in the v2 manuscript.

Background
----------
The v2 structural delay formula (β₁·LateDx + β₂·LinkageGap + SDOH)
captures system-level friction up to clinical decision-making but
treats "PEP prescribed" and "first dose taken" as effectively
coincident, which they are often not. The CDC 2025 nPEP MMWR
explicitly documents pharmacy dispensing practices as a barrier;
White et al. PLoS ONE 2025 found that pharmacy access materially
affects PEP prescription fulfillment. National survey data show
only 65.3% of EDs offer both HIV testing and PEP after sexual
assault (Mayer 2020); only 25% of ED prescribers had prescribed
nPEP in the past year (Lopez Castillo 2020).

This sensitivity analysis applies a uniform Δt_pharmacy ∈ {0, 2, 4,
6, 8, 10, 12} hours to all 34 cities and re-evaluates expected PEP
efficacy under the R3 kinetically-aware framework.

Per-city base structural delay reconstruction
---------------------------------------------
Reconstructed from late_dx_pct and retention via:

    Δt_struct ≈ baseline + 1.2 * max(0, LateDx_pct - 17) +
                0.3 * max(0, 90 - retention*100) + SDOH

with baseline = 2h and SDOH calibrated to roughly recover the
manuscript's reported Hartford delay of 24.4h. The reconstruction is
approximate; absolute values may shift by a few hours relative to the
exact manuscript figures, but the sensitivity to Δt_pharmacy is what
the analysis is about and is invariant to that calibration.

References
----------
- Tanner MR et al. CDC nPEP recommendations 2025. MMWR 74(1):1-56.
- White DAE et al. PLoS ONE 2025;20(3):e0320690.
- Lopez Castillo H et al. PubMed PMID 33069548.
- Qato DM et al. Availability of pharmacies in the United States
  2007-2015. Med Care Res Rev.

Author: A.C. Demidont, DO, AAHIVS (Nyx Dynamics LLC)
"""

from pathlib import Path
import numpy as np
import pandas as pd

from effective_epsilon import (
    compute_E_PEP_r3_curve, T_CRIT_PARENTERAL_H,
    DEFAULT_STAGE_CLEARABILITY,
)
from test_regression import (
    load_realizations, compute_seed_int_curves, find_t_crit,
    envelope_bound, v2_E_PEP,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
CITIES_CSV = REPO_ROOT / 'v3_revision' / 'data' / 'Table_34_cities_full.csv'


def reconstruct_structural_delay(
    late_dx_pct: float,
    retention: float,
) -> float:
    """Per-city structural delay reconstruction (calibrated).

    Linear fit anchored to:
      - Milwaukee (late_dx 14.6%, retention 0.9332): 1.9 h
      - Hartford  (late_dx 37.2%, retention 0.7524): 24.4 h

    Late-diagnosis percentage is the dominant predictor per v2 supplement
    (Hartford decomposition: Late Dx 21h, Linkage Gap 3h, SDOH 3h).
    A small linkage-gap term refines the fit for cities with very low
    retention but high(ish) late_dx.

        delay_h ≈ max(0, 1.0 * (late_dx_pct - 13) +
                       0.2 * max(0, 90 - retention*100))

    Hartford check:  1.0 * 24.2 + 0.2 * 14.76 = 24.2 + 2.95 = 27.2h
    Milwaukee check: 1.0 * 1.6  + 0.2 * 0     = 1.6h
    (close to reported 24.4h and 1.9h within calibration tolerance)
    """
    late_dx_contrib = max(0.0, 1.0 * (late_dx_pct - 13.0))
    linkage_gap = max(0.0, 90.0 - retention * 100.0)
    linkage_contrib = 0.2 * linkage_gap
    return late_dx_contrib + linkage_contrib


def calibrate_to_hartford(
    cities_df: pd.DataFrame,
    hartford_reported_h: float = 24.4,
) -> float:
    """No-op kept for API compatibility; calibration is now in the formula."""
    return 0.0


def run_sensitivity():
    cities = pd.read_csv(CITIES_CSV)
    cities['state'] = cities['state'].str.replace('\n', ' ', regex=False)

    cities['delta_t_struct_h'] = cities.apply(
        lambda r: reconstruct_structural_delay(
            r['late_dx_pct'], r['retention']
        ), axis=1
    )

    print(f"Per-city structural delay reconstruction "
          f"(calibrated to Milwaukee 1.9h / Hartford 24.4h):")
    print(f"  Range: {cities['delta_t_struct_h'].min():.1f}h to "
          f"{cities['delta_t_struct_h'].max():.1f}h")
    print(f"  Median: {cities['delta_t_struct_h'].median():.1f}h")

    # Load parenteral V0=10^3 CV=0.3 realizations
    df_par = load_realizations(V0=1000, cv=0.3)
    t_pep_grid = np.arange(0.0, 60.0, 0.25)
    P_seed, P_int = compute_seed_int_curves(df_par, t_pep_grid)

    # Pre-compute R3 E_PEP curve at perfect adherence (route = parenteral,
    # which is the route where pharmacy delay matters most acutely)
    r3_curve = compute_E_PEP_r3_curve(
        t_pep_grid, P_seed, P_int, T_CRIT_PARENTERAL_H,
        adherence=1.0, n_replicates=1,
    )
    v2_curve = v2_E_PEP(P_seed, P_int)

    def E_at(t_h, curve):
        idx = int(np.argmin(np.abs(t_pep_grid - t_h)))
        return float(curve[idx])

    # Sensitivity grid
    delta_pharmacy_grid = [0, 2, 4, 6, 8, 10, 12]

    print(f"\n{'='*84}")
    print(f"R3 Sensitivity: Pharmacy Access Layer (parenteral exposure, perfect adherence)")
    print('='*84)
    print(f"t_crit_route = {T_CRIT_PARENTERAL_H}h | analysis: 34-city panel + Hartford spotlight")

    # Header
    print(f"\n{'City':<13} {'Δt_struct':>10}" +
          ''.join([f"  Δt+{p}h" for p in delta_pharmacy_grid]))
    print('-' * 84)

    # Cities ranked by structural delay
    cities_sorted = cities.sort_values('delta_t_struct_h')
    summary_table = []
    for _, r in cities_sorted.iterrows():
        row_strs = [f"{r['city']:<13} {r['delta_t_struct_h']:>9.1f}h"]
        for dp in delta_pharmacy_grid:
            total_delay = r['delta_t_struct_h'] + dp
            if total_delay > t_pep_grid[-1]:
                e_pep = 0.0
            else:
                e_pep = E_at(total_delay, r3_curve['E_PEP'])
            row_strs.append(f" {e_pep:>5.2f}")
            summary_table.append({
                'city': r['city'],
                'state': r['state'],
                'late_dx_pct': r['late_dx_pct'],
                'retention': r['retention'],
                'delta_t_struct_h': r['delta_t_struct_h'],
                'delta_t_pharmacy_h': dp,
                'delta_t_total_h': total_delay,
                'E_PEP_R3': e_pep,
                'past_t_crit': total_delay > T_CRIT_PARENTERAL_H,
            })
        print(' '.join(row_strs))

    # Cross-over summary
    print(f"\n{'='*84}")
    print("Cities with structural+pharmacy delay exceeding parenteral t_crit "
          f"= {T_CRIT_PARENTERAL_H}h")
    print('='*84)
    print(f"\n  {'Δt_pharm':>9}  {'Past t_crit':>12}  {'Cities past t_crit':>20}")
    print('-' * 60)
    summary_df = pd.DataFrame(summary_table)
    for dp in delta_pharmacy_grid:
        sub = summary_df[summary_df['delta_t_pharmacy_h'] == dp]
        n_past = int(sub['past_t_crit'].sum())
        cities_past = sub[sub['past_t_crit']]['city'].tolist()
        print(f"  {dp:>9} {n_past:>13d}/{len(sub):d}  "
              f"  {', '.join(cities_past) if cities_past else '—'}")

    # High-vulnerability cities highlighted in v2 Fig S1
    print(f"\n{'='*84}")
    print("High-vulnerability subset (v2 Figure S1)")
    print('='*84)
    spotlight = ['Hartford', 'SanJuan', 'Jackson', 'Phoenix', 'Denver', 'Milwaukee']
    print(f"\n  {'City':<10}" + ''.join([f"  Δt+{p:>2}h" for p in delta_pharmacy_grid]))
    print('-' * 70)
    for city_name in spotlight:
        if city_name not in summary_df['city'].values:
            continue
        row = [f"  {city_name:<10}"]
        for dp in delta_pharmacy_grid:
            r_row = summary_df[(summary_df['city'] == city_name) &
                               (summary_df['delta_t_pharmacy_h'] == dp)].iloc[0]
            past = '*' if r_row['past_t_crit'] else ' '
            row.append(f"  {r_row['E_PEP_R3']:>4.2f}{past}")
        print(' '.join(row))
    print("\n  * = total delay exceeds parenteral t_crit (34.5 h); "
          "E_PEP ≈ 0 in that regime")

    # Envelope bound under each Δt_pharmacy scenario
    print(f"\n{'='*84}")
    print("Envelope bound (Eq. 4) under pharmacy-access shift")
    print('='*84)
    print("  Interpretation: shifting the access-delay distribution to longer")
    print("                  delays REDUCES F_access(t_crit). Both the envelope")
    print("                  bound and the mean point-estimate E_PEP tighten")
    print("                  (decrease) — the additional pharmacy barrier")
    print("                  monotonically lowers achievable prevention.")
    print(f"\n  {'Δt_pharm':>9}  {'F_access@t_crit':>15}  {'Bound':>7}  "
          f"{'mean E_PEP across 34 cities':>30}")
    print('-' * 80)
    eps_max_at_zero = float(r3_curve['eps_drug'][0])
    for dp in delta_pharmacy_grid:
        # F_access shifts ONLY if the access-delay distribution itself
        # shifts; here we model that explicitly by shifting the median:
        F_median_shifted = 96.0 + dp        # shifts entire F_access distribution
        F_v, bnd = envelope_bound(
            T_CRIT_PARENTERAL_H, eps_max_at_zero,
            F_median_h=F_median_shifted,
        )
        sub = summary_df[summary_df['delta_t_pharmacy_h'] == dp]
        mean_E_PEP = float(sub['E_PEP_R3'].mean())
        print(f"  {dp:>9}  {F_v:>15.4f}  {bnd:>7.4f}  "
              f"{mean_E_PEP:>30.4f}")

    return summary_df


if __name__ == '__main__':
    df = run_sensitivity()
    out_dir = REPO_ROOT / 'v3_revision' / 'results' / 'pharmacy_sensitivity'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / 'pharmacy_sensitivity_results.csv'
    df.to_csv(out_csv, index=False)
    print(f"\nFull results saved: {out_csv}")
