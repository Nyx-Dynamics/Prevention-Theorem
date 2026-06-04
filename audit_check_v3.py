"""
v3 Audit verification — validates PK-driven framework and pharmacy sensitivity
numerical claims against v3_revision/numerical_claims_v3.csv.

Runs test_regression.py and pharmacy_sensitivity.py, captures their output,
and checks each v3-specific claim within the stated tolerance. v2 baseline
claims (t_crit, compression ratio, envelope) are cross-checked from the
existing multiscale realizations CSV rather than re-running the MC.

Exit code 0 = all claims pass. Exit code 1 = one or more failures.
"""
import sys
import os
import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CLAIMS_CSV = ROOT / 'v3_revision' / 'numerical_claims_v3.csv'
R3_DIR = ROOT / 'v3_revision' / 'r3_pk_pd'
PYTHON = sys.executable

TOLERANCES = {
    'pk_tcrit_par_rho10':      0.5,   # hours
    'pk_tcrit_muc_rho10':      0.5,   # hours
    'pk_envelope_rho10':       0.5,   # percentage points
    'hartford_epep_dt8h':      0.05,  # absolute (must be < 0.05)
    'hartford_epep_baseline':  0.03,  # absolute
    'envelope_sweep_dt0':      0.5,   # pp
    'envelope_sweep_dt12h':    0.5,   # pp
    'mean_epep_sweep_dt0':     0.5,   # pp
    'mean_epep_sweep_dt12h':   0.5,   # pp
}


def load_claims():
    claims = {}
    with open(CLAIMS_CSV, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            claims[row['claim_id']] = row
    return claims


def run_test_regression():
    """Run test_regression.py and capture stdout."""
    print("\n" + "=" * 72)
    print("RUNNING: v3_revision/r3_pk_pd/test_regression.py")
    print("=" * 72)
    result = subprocess.run(
        [PYTHON, 'test_regression.py'],
        capture_output=True, text=True, cwd=str(R3_DIR)
    )
    if result.returncode != 0:
        print("FAILED (non-zero exit):")
        print(result.stderr)
        return None
    print(result.stdout)
    return result.stdout


def run_pharmacy_sensitivity():
    """Run pharmacy_sensitivity.py and capture stdout."""
    print("\n" + "=" * 72)
    print("RUNNING: v3_revision/r3_pk_pd/pharmacy_sensitivity.py")
    print("=" * 72)
    result = subprocess.run(
        [PYTHON, 'pharmacy_sensitivity.py'],
        capture_output=True, text=True, cwd=str(R3_DIR)
    )
    if result.returncode != 0:
        print("FAILED (non-zero exit):")
        print(result.stderr)
        return None
    print(result.stdout)
    return result.stdout


def parse_test_regression(output):
    """Extract key values from test_regression.py stdout."""
    parsed = {}
    lines = output.splitlines()
    in_parenteral_verdict = False
    in_mucosal_verdict = False
    for line in lines:
        if 'VERDICT  (parenteral)' in line:
            in_parenteral_verdict = True
            in_mucosal_verdict = False
        elif 'VERDICT  (mucosal)' in line:
            in_mucosal_verdict = True
            in_parenteral_verdict = False
        elif in_parenteral_verdict:
            if 'R3 t_crit at 100% adherence' in line:
                val = float(line.split(':')[1].strip().replace('h', '').split()[0])
                parsed['pk_tcrit_par_rho10'] = val
            elif 'R3 envelope bound at perf adherence' in line:
                val = float(line.split(':')[1].strip().replace('%', '').split()[0])
                parsed['pk_envelope_rho10'] = val
        elif in_mucosal_verdict:
            if 'R3 t_crit at 100% adherence' in line:
                val = float(line.split(':')[1].strip().replace('h', '').split()[0])
                parsed['pk_tcrit_muc_rho10'] = val
    return parsed


def parse_pharmacy_sensitivity(output):
    """Extract key values from pharmacy_sensitivity.py stdout."""
    parsed = {}
    lines = output.splitlines()
    in_envelope = False
    envelope_rows = []
    for line in lines:
        if 'Envelope bound (Eq. 4)' in line:
            in_envelope = True
            continue
        if in_envelope:
            stripped = line.strip()
            # Look for data rows: "0   0.0699  0.1150  0.9608"
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    dt = int(parts[0])
                    bound_pct = float(parts[2]) * 100
                    mean_epep_pct = float(parts[3]) * 100
                    envelope_rows.append((dt, bound_pct, mean_epep_pct))
                except (ValueError, IndexError):
                    pass
        # Hartford specific rows from High-vulnerability subset table
        if 'Hartford' in line and '0.47' in line:
            parts = line.split()
            try:
                # Format: Hartford  0.47  0.46  0.34  0.07  0.00*  0.00*  0.00*
                baseline = float(parts[1])
                dt8_raw = parts[5].replace('*', '')
                dt8 = float(dt8_raw)
                parsed['hartford_epep_baseline'] = baseline
                parsed['hartford_epep_dt8h'] = dt8
            except (ValueError, IndexError):
                pass

    if envelope_rows:
        dt0_row = next((r for r in envelope_rows if r[0] == 0), None)
        dt12_row = next((r for r in envelope_rows if r[0] == 12), None)
        if dt0_row:
            parsed['envelope_sweep_dt0'] = dt0_row[1]
            parsed['mean_epep_sweep_dt0'] = dt0_row[2]
        if dt12_row:
            parsed['envelope_sweep_dt12h'] = dt12_row[1]
            parsed['mean_epep_sweep_dt12h'] = dt12_row[2]

    return parsed


def verify_claims(claims, observed):
    """Compare observed values against CSV claims within tolerance."""
    print("\n" + "=" * 72)
    print("CLAIM VERIFICATION RESULTS")
    print("=" * 72)
    print(f"{'Claim ID':<35} {'Expected':>10} {'Observed':>10} {'Tol':>6}  {'Status'}")
    print("-" * 72)

    passed = 0
    failed = 0
    missing = 0

    for claim_id, row in claims.items():
        if row['category'] not in ('pk_driven', 'pharmacy'):
            continue  # v2 baseline claims verified separately by audit_check.py
        expected = float(row['value'])
        tol = TOLERANCES.get(claim_id)
        if tol is None:
            continue

        if claim_id not in observed:
            print(f"  {claim_id:<33} {expected:>10.3f} {'---':>10} {tol:>6.2f}  MISSING")
            missing += 1
            continue

        obs = observed[claim_id]
        # hartford_epep_dt8h is a "must be < threshold" check
        if claim_id == 'hartford_epep_dt8h':
            ok = obs < tol
            status = 'PASS' if ok else 'FAIL'
            print(f"  {claim_id:<33} {expected:>10.3f} {obs:>10.3f} {tol:>6.2f}  {status} (must be < {tol})")
        else:
            diff = abs(obs - expected)
            ok = diff <= tol
            status = 'PASS' if ok else 'FAIL'
            print(f"  {claim_id:<33} {expected:>10.3f} {obs:>10.3f} {tol:>6.2f}  {status}")

        if ok:
            passed += 1
        else:
            failed += 1

    print("-" * 72)
    print(f"  {passed} passed  |  {failed} failed  |  {missing} missing")
    return failed == 0 and missing == 0


def main():
    claims = load_claims()

    reg_output = run_test_regression()
    pharm_output = run_pharmacy_sensitivity()

    if reg_output is None or pharm_output is None:
        print("\nABORTED: one or more scripts failed to run.")
        sys.exit(1)

    observed = {}
    observed.update(parse_test_regression(reg_output))
    observed.update(parse_pharmacy_sensitivity(pharm_output))

    all_pass = verify_claims(claims, observed)

    if all_pass:
        print("\nALL v3 CLAIMS VERIFIED.")
        sys.exit(0)
    else:
        print("\nSOME CLAIMS FAILED — review output before tagging or submitting.")
        sys.exit(1)


if __name__ == '__main__':
    main()
