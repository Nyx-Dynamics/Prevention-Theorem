"""
Regression test for the envelope-bound correction.

Two things must hold simultaneously:

  1. The PREVIOUSLY REPORTED single value (11.3% at the canonical
     operating point) is still reproducible, unchanged, via
     envelope_bound() — the manuscript-reproduction pathway is intact.

  2. The value reported in the CROI abstract (6.6-11.3%) is produced by
     envelope_bound_range(), whose upper end is identical to (1).

The correction is additive: no existing function changed behaviour, and
no t_crit or Monte Carlo output is touched by this module.

Run:  python3 test_envelope_bound_range.py
"""
import sys

from test_regression import (
    envelope_bound, envelope_bound_range, envelope_band,
)

# Canonical operating point (manuscript Eq. 4 / numerical_claims_v4.csv)
W_CANONICAL = 34.5      # h, parenteral t_crit at eta=0.05
EPS_MAX = 0.95
F_MEDIAN = 96.0         # h, log-normal access-delay median
GSD = 2.0
ETA = 0.05

TOL = 5e-4              # 0.05 percentage points

failures = []


def check(label, got, want, tol=TOL):
    ok = abs(got - want) <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {label:<58} got {got*100:7.3f}%  want {want*100:7.3f}%")
    if not ok:
        failures.append(label)
    return ok


print("=" * 92)
print("ENVELOPE BOUND — additive correction regression")
print("=" * 92)

# ---------------------------------------------------------------- 1. legacy
print("\n1. Legacy single-value pathway is unchanged (upper-bound convention)")
F_legacy, bound_legacy = envelope_bound(W_CANONICAL, EPS_MAX, F_MEDIAN, GSD, ETA)
check("envelope_bound() F_access(34.5h)", F_legacy, 0.0699, tol=5e-4)
check("envelope_bound() reproduces the reported 11.3%", bound_legacy, 0.113)

# ---------------------------------------------------------------- 2. range
print("\n2. CROI-reported range comes from envelope_bound_range()")
F_rng, lower, upper = envelope_bound_range(W_CANONICAL, EPS_MAX, F_MEDIAN, GSD, ETA)
check("envelope_bound_range() F_access matches legacy", F_rng, F_legacy, tol=1e-12)
check("lower end (no residual protection)", lower, 0.066)
check("upper end (residual efficacy up to eta)", upper, 0.113)

# ---------------------------------------------------------------- 3. identity
print("\n3. Upper end is IDENTICAL to the legacy value, not merely close")
same = (upper == bound_legacy)
print(f"  {'PASS' if same else 'FAIL'}  upper == envelope_bound()"
      f"{'':<34} {upper!r} vs {bound_legacy!r}")
if not same:
    failures.append("upper == legacy")

# ---------------------------------------------------------------- 4. ordering
print("\n4. Structural invariants")
inv = [("lower <= upper", lower <= upper),
       ("lower == F_access * eps_max", abs(lower - F_rng * EPS_MAX) < 1e-12),
       ("upper - lower == (1-F)*eta", abs((upper - lower) - (1 - F_rng) * ETA) < 1e-12),
       ("eta=0 collapses range", envelope_bound_range(W_CANONICAL, EPS_MAX, F_MEDIAN, GSD, 0.0)[1]
        == envelope_bound_range(W_CANONICAL, EPS_MAX, F_MEDIAN, GSD, 0.0)[2])]
for label, ok in inv:
    print(f"  {'PASS' if ok else 'FAIL'}  {label}")
    if not ok:
        failures.append(label)

# ---------------------------------------------------------------- 5. CROI band
print("\n5. CROI graphic band (generated from envelope_bound_range)")
band = envelope_band([28, 34.5, 48, 60, 72], EPS_MAX, F_MEDIAN, GSD, ETA)
print(f"\n  {'W (h)':>7}{'F_access':>11}{'lower':>9}{'upper':>9}   published in CROI table")
EXPECTED = {28: (0.038, 3.6, 8.4), 34.5: (0.070, 6.6, 11.3), 48: (0.159, 15.1, 19.3),
            60: (0.249, 23.6, 27.4), 72: (0.339, 32.2, 35.5)}
for _, r in band.iterrows():
    W = r['window_h']
    eF, elo, ehi = EXPECTED[W]
    ok = (abs(r['F_access'] - eF) < 5e-4 and abs(r['ceiling_lower_pct'] - elo) < 0.05
          and abs(r['ceiling_upper_pct'] - ehi) < 0.05)
    print(f"  {W:>7}{r['F_access']:>11.4f}{r['ceiling_lower_pct']:>8.1f}%"
          f"{r['ceiling_upper_pct']:>8.1f}%   {elo}-{ehi}%  {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(f"band W={W}")

print()
print("=" * 92)
if failures:
    print(f"FAILED ({len(failures)}): " + "; ".join(failures))
    sys.exit(1)
print("ALL PASS — legacy 11.3% reproduces unchanged; CROI range 6.6-11.3% "
      "produced by envelope_bound_range()")
sys.exit(0)
