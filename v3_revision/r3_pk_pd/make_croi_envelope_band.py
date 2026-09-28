"""
Generate the CROI 2027 graphic artifact: the window-sensitivity band of the
distribution-free F_access x t_crit bound.

Produces v3_revision/results/croi_2027/croi_envelope_band.csv, the artifact
cited by numerical_claims_v4.csv entry `croi_envelope_band_5row`.

Both ends of the band come from envelope_bound_range():

    lower = F_access(W) * eps_max
              zero residual post-window protection (v4 gating rule: no drug
              acquired before t_crit implies E_PEP = 0)

    upper = lower + (1 - F_access(W)) * eta
              preserves the prior eta = 0.05 residual-efficacy convention;
              identical to the legacy envelope_bound()

This script writes an artifact only.  It does not touch the
manuscript-reproduction pathway, t_crit, or any Monte Carlo output.

Run:  python3 make_croi_envelope_band.py
"""
from pathlib import Path

from test_regression import envelope_band, envelope_bound, envelope_bound_range

# Canonical parameterization (manuscript Eq. 4 / numerical_claims_v4.csv)
WINDOWS_H = [28, 34.5, 48, 60, 72]
EPS_MAX = 0.95
F_MEDIAN_H = 96.0
GSD = 2.0
ETA = 0.05
W_CANONICAL = 34.5

OUT_DIR = Path(__file__).resolve().parents[1] / 'results' / 'croi_2027'
OUT_CSV = OUT_DIR / 'croi_envelope_band.csv'


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    band = envelope_band(WINDOWS_H, EPS_MAX, F_MEDIAN_H, GSD, ETA)
    band['eps_max'] = EPS_MAX
    band['F_median_h'] = F_MEDIAN_H
    band['GSD'] = GSD
    band['eta'] = ETA
    band['is_canonical_window'] = band['window_h'] == W_CANONICAL
    band['lower_convention'] = 'zero residual post-window protection'
    band['upper_convention'] = 'eta=0.05 residual efficacy (legacy convention)'
    band.to_csv(OUT_CSV, index=False)

    print('=' * 88)
    print('CROI 2027 envelope band')
    print(f'  envelope_band(windows_h={WINDOWS_H}, eps_max_at_zero={EPS_MAX}, '
          f'F_median_h={F_MEDIAN_H}, GSD={GSD}, eta={ETA})')
    print('=' * 88)
    print(f"  {'W (h)':>7}{'F_access':>11}{'lower %':>10}{'upper %':>10}   convention")
    for _, r in band.iterrows():
        mark = '  <- canonical' if r['is_canonical_window'] else ''
        print(f"  {r['window_h']:>7}{r['F_access']:>11.4f}"
              f"{r['ceiling_lower_pct']:>10.1f}{r['ceiling_upper_pct']:>10.1f}"
              f"   lower = F*eps_max; upper = +(1-F)*eta{mark}")

    # Cross-check the canonical row against both callable pathways
    _, legacy = envelope_bound(W_CANONICAL, EPS_MAX, F_MEDIAN_H, GSD, ETA)
    _, lo, hi = envelope_bound_range(W_CANONICAL, EPS_MAX, F_MEDIAN_H, GSD, ETA)
    assert hi == legacy, 'upper end must equal the legacy envelope_bound()'
    row = band[band['window_h'] == W_CANONICAL].iloc[0]
    assert abs(row['ceiling_lower_pct'] - lo * 100) < 1e-9
    assert abs(row['ceiling_upper_pct'] - hi * 100) < 1e-9

    print()
    print(f'  canonical W={W_CANONICAL}h -> {lo*100:.1f}% - {hi*100:.1f}% '
          f'(upper is bit-identical to legacy envelope_bound(): {hi == legacy})')
    print(f'  saved {OUT_CSV}')


if __name__ == '__main__':
    main()
