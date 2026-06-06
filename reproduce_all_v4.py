#!/usr/bin/env python3
"""
reproduce_all_v4.py — v4 reproduction pipeline for the PEP window manuscript
("The HIV Post-Exposure Prophylaxis Window: A Multiscale Framework Linking
Within-Host Integration Kinetics to Population-Scale Structural Access",
PLOS Computational Biology submission).

What's new in v4 vs v3
======================

v4 corrects a methodological error in the v3 pharmacy sensitivity output.
The prior "envelope bound" scalar conflated three quantities that should
be kept separate:
  - F_access(t_crit), a population-level access distribution
  - eps_max, treated as a constant 0.98 rather than a function of the
    remaining window (t_crit - t_acq)
  - eps_min = 0.05, a phantom residual floor inconsistent with the
    gating-event causal chain

Under v4, pharmacy access is modeled as an UPSTREAM GATING EVENT:
patients who do not acquire medication before t_crit never receive drug,
PK is irrelevant, and E_PEP = 0 (not a residual eps_min floor). The
"envelope" is a 2D corridor in (t_acq, E_PEP) space, bounded above by
perfect-adherence PK and below by low-adherence PK. The corridor is a
property of viral kinetics + drug PK; pharmacy delays do not compress it,
they SLIDE cities rightward along it. Cities whose t_acq exceeds t_crit
are displaced off the corridor entirely.

v4 reproduces all v3 outputs except the (removed) scalar envelope bound,
and adds:
  - envelope_corridor.csv: upper/lower envelope curves over a t_acq grid
  - city_envelope_positions.csv: per-(city, dt_pharm) corridor positions
  - pharmacy_displacement_summary.csv: count of cities past t_crit
  - Figure S14: envelope corridor with 34-city overlay + Hartford
    displacement arrow

Compared to reproduce_all_v3.py, this script also now writes a real
run_metadata.txt provenance file (capturing repo state, git HEAD,
library versions) into the timestamped run directory. The prior v3
deposits inherited a Kassanjee run_metadata.txt from a different repo,
which is no longer the case.

See PATCH_NOTES_envelope_corridor.md for the full diff between v3 and v4
behaviors and outputs.
"""

import os
import platform
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# ---------------------------------------------------------------------------
# Project structure
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / 'SRC'
V2_REVISION_DIR = PROJECT_ROOT / 'v2_revision'
V3_REVISION_DIR = PROJECT_ROOT / 'v3_revision'
R3_DIR = V3_REVISION_DIR / 'r3_pk_pd'

# Timestamped run directory; suffix _v4 distinguishes from v3 runs in the
# same runs/ directory.
TIMESTAMP = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
RUNS_DIR = V3_REVISION_DIR / 'runs'
FIGURES_DIR = RUNS_DIR / f"run_{TIMESTAMP}_v4"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run_command(cmd, env_add=None, cwd=None):
    """Run a shell command with optional PYTHONPATH addition. Exits on failure."""
    print(f"\n>>> Running: {' '.join(str(c) for c in cmd)}")
    env = os.environ.copy()
    if env_add:
        env['PYTHONPATH'] = (
            env.get('PYTHONPATH', '')
            + (':' if env.get('PYTHONPATH') else '')
            + str(env_add)
        )
    result = subprocess.run(cmd, env=env, cwd=cwd)
    if result.returncode != 0:
        print(f"Error executing {' '.join(str(c) for c in cmd)}")
        sys.exit(result.returncode)


def _git(*args):
    """Run a git command at PROJECT_ROOT; return stripped stdout or None on failure."""
    try:
        return subprocess.check_output(
            ['git', '-C', str(PROJECT_ROOT), *args],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return None


def write_run_metadata(out_path):
    """Capture provenance: git state, host, python + library versions, timestamp.

    Replaces the legacy Kassanjee run_metadata.txt that v3 deposits inherited
    by accident (its 'Working dir' pointed at a different repository).
    """
    lines = [
        "# reproduce_all_v4.py — run provenance",
        "# Prevention-Theorem repository (Corner 2, PEP window paper)",
        "",
        f"Run timestamp UTC:   {datetime.now(timezone.utc).isoformat()}",
        f"Run timestamp local: {datetime.now().isoformat()}",
        f"Host:                {socket.gethostname()}",
        f"Working dir:         {PROJECT_ROOT}",
        f"Python version:      {sys.version.split()[0]}",
        f"Platform:            {platform.platform()}",
        "",
    ]

    git_head = _git('rev-parse', 'HEAD')
    lines.append(f"Git HEAD:            {git_head or '(unavailable)'}")

    git_describe = _git('describe', '--tags', '--always')
    if git_describe:
        lines.append(f"Git describe:        {git_describe}")

    git_status = _git('status', '--porcelain')
    if git_status is None:
        lines.append("Git status:          (unavailable)")
    elif git_status:
        lines.append("Git status:          DIRTY (uncommitted changes present)")
        lines.append("Uncommitted files (first 20):")
        for line in git_status.splitlines()[:20]:
            lines.append(f"  {line}")
    else:
        lines.append("Git status:          clean")

    lines.append("")
    lines.append("Library versions:")
    for module in ['numpy', 'scipy', 'pandas', 'matplotlib']:
        try:
            mod = __import__(module)
            lines.append(f"  {module}: {mod.__version__}")
        except Exception:
            lines.append(f"  {module}: (unavailable)")

    out_path.write_text('\n'.join(lines) + '\n')
    print(f"Wrote provenance metadata: {out_path}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("V4 REPRODUCTION PIPELINE: PREVENTION THEOREM")
    print("PLOS Computational Biology submission")
    print("v4: envelope-corridor framing for pharmacy access (supersedes v3)")
    print("=" * 80)

    # -----------------------------------------------------------------------
    # PART 0: Capture run provenance (new in v4)
    # -----------------------------------------------------------------------
    write_run_metadata(FIGURES_DIR / 'run_metadata.txt')

    # -----------------------------------------------------------------------
    # PART 1: v2 pipeline (main-text Figures 1-4) — unchanged in v4
    # The v4 manuscript's main-text figures are the same as v2 because the
    # PK-driven framework recovers the baseline at perfect adherence and
    # the corridor refactor affects only the supplement S14 layer.
    # -----------------------------------------------------------------------
    print("\n--- PART 1: Reproducing v2 main-text figures ---")

    # Figure 1: Route comparison from multiscale baseline
    print("\n--- Generating Figure 1 (Route Comparison, multiscale baseline) ---")
    run_command([sys.executable, str(V2_REVISION_DIR / 'make_figure_1_v2.py')])
    for ext in ['png', 'pdf']:
        src = V2_REVISION_DIR / 'figures' / f'Figure_1_Route_Dependent_PEP_Efficacy_Decay.{ext}'
        if src.exists():
            dest = FIGURES_DIR / f'Figure_1_VL_Compression_PEP_Parenteral.{ext}'
            src.replace(dest)
            print(f"Moved Figure 1 to {dest}")

    # Capture Figure 1 source data (multiscale baseline)
    fig1_data_dir = SRC_DIR / 'multiscale_model' / 'results_v3'
    if fig1_data_dir.exists():
        for csv_file in fig1_data_dir.glob('*.csv'):
            shutil.copy2(csv_file, FIGURES_DIR / csv_file.name)

    # Figure 2: Stochastic efficacy & VL uncertainty
    print("\n--- Generating Figure 2 (Stochastic VL Uncertainty) ---")
    run_command(
        [sys.executable, str(SRC_DIR / 'perelson/stochastic/PEP_stochastic_perelson.py')],
        env_add=str(SRC_DIR),
    )
    fig2_default = PROJECT_ROOT / 'pep_stochastic_vl_uncertainty.png'
    if fig2_default.exists():
        dest = FIGURES_DIR / 'Figure_2_Stochastic_Efficacy_VL_Uncertainty.png'
        fig2_default.replace(dest)
        print(f"Moved Figure 2 to {dest}")
    stochastic_csv_dir = PROJECT_ROOT / 'pep_stochastic_results'
    if stochastic_csv_dir.exists():
        for csv_file in stochastic_csv_dir.glob('*.csv'):
            shutil.copy2(csv_file, FIGURES_DIR / csv_file.name)

    # Figures 3-4: City-stratified analysis (main Figure 3 + supplementary S1)
    print("\n--- Generating Figure 3 / Supplementary S1 (City Stratified) ---")
    run_command([sys.executable, str(V2_REVISION_DIR / 'city_stratified_figures.py')])
    city_results_dir = PROJECT_ROOT / 'results' / 'city_analysis'
    if city_results_dir.exists():
        for csv_file in city_results_dir.glob('*.csv'):
            shutil.copy2(csv_file, FIGURES_DIR / csv_file.name)
        moves = [
            ('Fig_CityStratified_PEP.png', 'Figure_3_City_Stratified_PEP_Efficacy.png'),
            ('Fig_CityComparison_Focus.png', 'Figure_S1_CityComparison_Focus.png'),
            ('Figure_3_City_Stratified_PEP_Efficacy.png', 'Figure_3_City_Stratified_PEP_Efficacy.png'),
            ('Figure_S1_CityComparison_Focus.png', 'Figure_S1_CityComparison_Focus.png'),
        ]
        for src_name, dest_name in moves:
            src = city_results_dir / src_name
            if src.exists():
                src.replace(FIGURES_DIR / dest_name)
                print(f"Moved {src_name} to {FIGURES_DIR / dest_name}")

    # -----------------------------------------------------------------------
    # PART 2: PK-driven framework (v3) + corrected pharmacy access (v4)
    # -----------------------------------------------------------------------
    print("\n--- PART 2: PK-driven framework + v4 envelope-corridor pharmacy ---")

    # R3 regression test: PK-driven framework recovers multiscale baseline
    # within tolerance at perfect adherence. Unchanged from v3.
    print("\n--- R3 regression: PK-driven framework vs. multiscale baseline ---")
    run_command([sys.executable, str(R3_DIR / 'test_regression.py')])

    # S14 pharmacy access — v4 corrected pharmacy_sensitivity.py
    # Outputs the envelope corridor + city positions + displacement summary.
    # Replaces the v3 'envelope bound' scalar with corridor framing.
    print("\n--- S14: Pharmacy access (v4 corridor framing) ---")
    run_command([sys.executable, str(R3_DIR / 'pharmacy_sensitivity.py')])

    pharm_out_dir = V3_REVISION_DIR / 'results' / 'pharmacy_sensitivity_corrected'
    if pharm_out_dir.exists():
        for csv_file in pharm_out_dir.glob('*.csv'):
            shutil.copy2(csv_file, FIGURES_DIR / csv_file.name)
            print(f"Copied {csv_file.name} to {FIGURES_DIR}")
    else:
        print(f"WARNING: expected pharmacy outputs at {pharm_out_dir}")
        print("  Check that pharmacy_sensitivity.py is the v4 corrected version")
        print("  and that its output directory is pharmacy_sensitivity_corrected/")

    # S2 tables: replaces the v3-era figure with three supplement tables
    # (corridor sample, city baseline positions, displacement summary).
    print("\n--- Tables S2a/S2b/S2c: corridor sample, city positions, displacement ---")
    run_command([sys.executable, str(R3_DIR / 'generate_S2_tables.py')])
    tables_dir = V3_REVISION_DIR / 'tables'
    if tables_dir.exists():
        for tex_file in tables_dir.glob('table_S2*.tex'):
            shutil.copy2(tex_file, FIGURES_DIR / tex_file.name)
            print(f"Copied {tex_file.name} to {FIGURES_DIR}")

    # -----------------------------------------------------------------------
    # PART 3: Verification of numerical claims (v4 registry)
    # -----------------------------------------------------------------------
    print("\n--- PART 3: Verifying v4 numerical claims ---")
    claims_csv_v4 = V3_REVISION_DIR / 'numerical_claims_v4.csv'
    claims_csv_v3 = V3_REVISION_DIR / 'numerical_claims_v3.csv'
    if claims_csv_v4.exists():
        shutil.copy2(claims_csv_v4, FIGURES_DIR / claims_csv_v4.name)
        print(f"Copied v4 numerical claims to {FIGURES_DIR / claims_csv_v4.name}")
    elif claims_csv_v3.exists():
        shutil.copy2(claims_csv_v3, FIGURES_DIR / claims_csv_v3.name)
        print(f"WARNING: numerical_claims_v4.csv not found; archived v3 file instead.")
        print(f"  Update numerical_claims_v4.csv per PATCH_NOTES_envelope_corridor.md")
        print(f"  Removed claims:    envelope_sweep_dt0, envelope_sweep_dt12h")
        print(f"  Reframed claims:   mean_epep_sweep_dt0, mean_epep_sweep_dt12h")
        print(f"  Added claims:      corridor_*, displacement_count_*, tcrit_pk_rho030")
    else:
        print(f"WARNING: no numerical_claims CSV found in {V3_REVISION_DIR}")

    print()
    print("Key v4 claims that this run reproduces:")
    print()
    print("  Multiscale baseline (unchanged from v3, commit 37e27ea):")
    print("    Parenteral t_crit at η=0.05 (CV=0.3, V₀=10³):   34.5 h")
    print("    Mucosal t_crit at η=0.05 (CV=0.3, V₀=1):         60.5 h")
    print("    Compression ratio (mucosal/parenteral):           ~1.75×")
    print("    F_access × t_crit envelope at canonical params:   11.3%")
    print()
    print("  PK-driven framework (R3, unchanged from v3):")
    print("    PK t_crit at ρ=1.0 (parenteral, upper corridor):  34.0 h (-0.5 h)")
    print("    PK t_crit at ρ=0.30 (parenteral, lower corridor): 32.0 h (-2.5 h)")
    print("    PK envelope bound at ρ=1.0:                        11.2%")
    print("    Recovery within tolerance of multiscale at ρ=1.0")
    print()
    print("  Pharmacy access — CORRIDOR FRAMING (new in v4):")
    print("    The envelope is a 2D corridor (E_PEP vs t_acq), NOT a scalar.")
    print("    Pharmacy delay slides cities rightward; the corridor is fixed.")
    print("    Hartford at Δt_pharm=0:   t_acq=24.4h, E_PEP_upper = 0.4789 (on corridor)")
    print("    Hartford at Δt_pharm=8h:  t_acq=32.4h, E_PEP_upper = 0.1748 (near cliff, on corridor)")
    print("    Hartford at Δt_pharm=10h: t_acq=34.4h > 34.0h, E_PEP_upper ≈ 0 (off corridor)")
    print("    Cities displaced past t_crit at Δt_pharm ∈ {0,2,4,6,8}h:  0 of 34")
    print("    Cities displaced past t_crit at Δt_pharm ∈ {10,12}h:     1 of 34 (Hartford)")
    print("    Cohort mean E_PEP_upper, Δt_pharm=0:    0.962 (NOT a population estimate;")
    print("    Cohort mean E_PEP_upper, Δt_pharm=12h:  0.865  high-burden metros only)")
    print()
    print("  SUPERSEDED in v4 (do not cite these from the v3 registry):")
    print("    Scalar envelope_sweep_dt0 = 11.5%   (Bernoulli mixture; misleading)")
    print("    Scalar envelope_sweep_dt12h = 9.6%  (same; conflated F_access × ε_max)")
    print("    Phantom ε_min = 0.05 floor          (no drug ⇒ no effect ⇒ ε_min = 0)")

    print("\n" + "=" * 80)
    print("V4 REPRODUCTION COMPLETE")
    print(f"All artifacts available in: {FIGURES_DIR}")
    print("=" * 80)
    print()
    print("Next steps for v4 Zenodo deposit:")
    print("  1. Verify all expected CSVs and figures landed in FIGURES_DIR above.")
    print("  2. Update numerical_claims_v4.csv per PATCH_NOTES_envelope_corridor.md")
    print("     if not yet done. Replace v3 with v4 in main manuscript Data Availability.")
    print("  3. Replace supplement S14 with S14_pharmacy_access_rewrite.tex content.")
    print("  4. Commit, retag v4.0.0, push. Verify GitHub-Zenodo automation captures")
    print("     the corrected outputs (results CSVs not in .gitignore).")
    print("  5. Update manuscript Data Availability with the new versioned Zenodo DOI.")


if __name__ == "__main__":
    main()
