"""
Audit verification — tests BOTH copies of the parenteral model
(project root and SRC/) to prove they produce identical t_crit values.

The Makefile runs scripts from SRC/, so SRC/PEP_parenteral_perelson.py is
the canonical version that produced the published figures. We test both
copies to confirm any divergence is cosmetic only.

Manuscript claims (parenteral, η=0.05): t_crit ≈ 16-28 h
Audit prediction (parenteral, η=0.05):  t_crit ≈ 60-68 h at PWID-relevant VL
The 16-28 h band is the seeding_midpoint, not t_crit.
"""
import sys, os, importlib
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC  = os.path.join(ROOT, 'SRC')

vl_grid = [50, 200, 1000, 10000, 30000, 100000, 1000000]
hours = np.linspace(0, 200, 2000)


def test_version(version_label, path_to_add):
    """Import ParenteralExposureModel from a specific path and compute
    t_crit values across the VL grid."""
    print()
    print("=" * 78)
    print(f"TESTING: {version_label}")
    print(f"  sys.path[0] = {path_to_add}")
    print("=" * 78)

    # Force re-import from this path
    for mod in ['PEP_parenteral_perelson']:
        if mod in sys.modules:
            del sys.modules[mod]
    sys.path.insert(0, path_to_add)

    from PEP_parenteral_perelson import ParenteralExposureModel
    print(f"  Loaded from: {sys.modules['PEP_parenteral_perelson'].__file__}")
    print()

    print(f"{'Source VL':<14} {'seeding_mid':<14} {'int_complete':<14} "
          f"{'t_crit η=0.05':<16} {'t_crit η=0.50':<16}")
    print(f"{'(copies/mL)':<14} {'(hours)':<14} {'(hours)':<14} "
          f"{'(hours)':<16} {'(hours)':<16}")
    print("-" * 76)

    rows = []
    for vl in vl_grid:
        m = ParenteralExposureModel(source_viral_load=vl,
                                     exposure_type='pwid_shared_needle')
        eff = np.array([m.pep_efficacy(h)['pep_efficacy'] for h in hours])

        below_05 = np.where(eff < 0.05)[0]
        below_50 = np.where(eff < 0.50)[0]
        t_05 = float(hours[below_05[0]]) if len(below_05) > 0 else float('nan')
        t_50 = float(hours[below_50[0]]) if len(below_50) > 0 else float('nan')

        print(f"{vl:<14} {m.seeding_midpoint:<14.1f} "
              f"{m.integration_complete:<14.1f} "
              f"{t_05:<16.1f} {t_50:<16.1f}")
        rows.append((vl, m.seeding_midpoint, m.integration_complete, t_05, t_50))

    # Clean up sys.path
    sys.path.remove(path_to_add)
    return rows


# Run both versions
print("=" * 78)
print("AUDIT VERIFICATION — Prevention-Theorem PEP model")
print("=" * 78)

rows_root = test_version("Project-root version (PEP_parenteral_perelson.py)", ROOT)
rows_src  = test_version("SRC/ version (Makefile uses this one)", SRC)

# Cross-check: do the two versions produce identical numbers?
print()
print("=" * 78)
print("CROSS-CHECK: project-root output == SRC/ output?")
print("=" * 78)
all_match = True
for r, s in zip(rows_root, rows_src):
    if r != s:
        all_match = False
        print(f"  MISMATCH at VL={r[0]}: root={r[1:]}, src={s[1:]}")
if all_match:
    print("  ✓ All 7 VL points produce identical t_crit values across both versions.")
    print("    The cosmetic differences in plotting code do not affect the model.")
else:
    print("  ✗ Versions diverge in numerical output. Investigate before continuing.")

# Mucosal verification
print()
print("=" * 78)
print("MUCOSAL VERIFICATION (SRC/PEP_mucosal.py)")
print("=" * 78)
sys.path.insert(0, SRC)
from PEP_mucosal import InfectionEstablishmentModel

m = InfectionEstablishmentModel()
print(f"  seeding_midpoint     = {m.seeding_midpoint} h  (line 25 of PEP_mucosal.py)")
print(f"  integration_complete = {m.integration_complete} h")
print()

eff = np.array([m.pep_efficacy(h)['pep_efficacy'] for h in hours])
print(f"  {'Threshold η':<15} {'t_crit (h)':<15}")
print(f"  {'-'*30}")
for eta in [0.80, 0.50, 0.10, 0.05, 0.01]:
    below = np.where(eff < eta)[0]
    t_crit = float(hours[below[0]]) if len(below) > 0 else float('nan')
    flag = "  ← manuscript claims 68-76h here" if eta == 0.05 else ""
    print(f"  η = {eta:<11.2f} {t_crit:<15.1f}{flag}")

print()
print("=" * 78)
print("INTERPRETATION")
print("=" * 78)
print("""
The manuscript text states:
    parenteral t_crit at η=0.05 ≈ 16-28 hours
    mucosal   t_crit at η=0.05 ≈ 68-76 hours

The audit predicts:
    parenteral t_crit at η=0.05 = 50-70 hours at PWID-relevant VL
    mucosal   t_crit at η=0.05 = ~133 hours
    The "16-28h" band corresponds to the parenteral SEEDING_MIDPOINT
      (the 'seeding_mid' column in the parenteral table).
    The "68-76h" band corresponds to the mucosal SEEDING_MIDPOINT (72h)
      or to mucosal t_crit at η=0.80 (about 66h).

If your output above matches the audit prediction, the manuscript text
conflates seeding_midpoint with t_crit at η=0.05. This is a real
discrepancy that needs to be reconciled before/during the revision round.

If your output disagrees with the audit (e.g., parenteral t_crit at
η=0.05 actually shows 16-28h on your machine), then something in your
local environment overrides the model behavior, and we need to find what.
""")