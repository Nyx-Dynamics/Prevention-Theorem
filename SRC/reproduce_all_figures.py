"""
reproduce_all_figures.py
========================
Master wrapper for generating all manuscript and supplement figures
and CSVs for the Science submission of the Prevention Theorem paper.

WHAT THIS PRODUCES
------------------
Figures:
  Fig 1  → pep_parenteral_vl_sweep.png          (Script 1: PEP_parenteral_perelson.py)
  Fig 2  → pep_stochastic_vl_uncertainty.png     (Script 2: PEP_stochastic_perelson.py)
  Fig 3  → Fig_CityStratified_PEP.png            (Script 4: city_stratified_figures.py)
  S3     → Fig_CityComparison_Focus.png          (Script 4: city_stratified_figures.py)

CSVs (stochastic results):
  pep_stochastic_efficacy_curves.csv
  pep_stochastic_summary_timepoints.csv
  pep_vl_knowledge_premium.csv
  pep_stochastic_ci_width_regimes.csv

CSVs (city analysis):
  city_vl_profiles.csv
  city_pep_efficacy_results.csv
  city_pep_efficacy_24h_counterfactual.csv
  city_structural_delay_cost.csv

KNOWN BUGS FIXED BY THIS WRAPPER (do not edit source scripts):
  1. city_stratified_figures.py lines 33-36 — hardcoded stale upload
     paths. This wrapper patches them to use freshly generated CSVs.
  2. city_stratified_figures.py lines 105-106 — legend labels wrong.
     Source says 'Moderate barrier (8-30h)' and 'High barrier (>30h)'.
     Correct per aidsvu_city_profiles.py thresholds (<8h, <18h, <30h, >=30h):
       Moderate barrier: 8-18h
       High barrier:     18-30h
  3. Makefile assumes SRC/ subdir — does not exist.
     All scripts run from project root.
  4. aidsvu_city_profiles.py writes to /mnt/user-data/outputs/ when run
     standalone. This wrapper captures the DataFrames directly in memory
     and passes them to the plotting script, avoiding path dependency.

USAGE
-----
  cd /path/to/Prevention-Theorem
  python reproduce_all_figures.py

  All outputs written to ./outputs/ (created if absent).
  Copy outputs/ contents to your submission package.

DEPENDENCIES
------------
  Python 3.8+
  numpy, scipy, matplotlib, pandas, openpyxl
  All *.py scripts and *_AIDSVu_*.xlsx files must be in the same directory
  as this script (the project root).

Zenodo: https://doi.org/10.5281/zenodo.18746065
GitHub: https://github.com/Nyx-Dynamics/Prevention-Theorem
Author: AC Demidont, DO, AAHIVS | Nyx Dynamics LLC | March 2026
"""

import os
import sys
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# ── resolve project root (same directory as this script) ─────────────────────
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'outputs')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── stochastic CSV subdirectory (matches PEP_stochastic_perelson.py default) ──
STOCHASTIC_CSV_DIR = os.path.join(PROJECT_ROOT, 'pep_stochastic_results')
os.makedirs(STOCHASTIC_CSV_DIR, exist_ok=True)

print("=" * 70)
print("PREVENTION THEOREM — REPRODUCE ALL FIGURES")
print(f"Project root : {PROJECT_ROOT}")
print(f"Output dir   : {OUTPUT_DIR}")
print("=" * 70)


# =============================================================================
# STEP 1 — Fig 1: PEP_parenteral_perelson.py
# Produces: pep_parenteral_vl_sweep.png
# No bugs in this script. Run it directly with the save path redirected.
# =============================================================================

def run_fig1():
    t0 = time.time()
    print("\n[1/4] Generating Fig 1 — Parenteral VL sweep…")

    from PEP_parenteral_perelson import plot_vl_sweep

    save_path = os.path.join(OUTPUT_DIR, 'pep_parenteral_vl_sweep.png')
    plot_vl_sweep(save_path=save_path)
    plt.close('all')

    print(f"      → {save_path}  ({time.time()-t0:.0f}s)")


# =============================================================================
# STEP 2 — Fig 2 + stochastic CSVs: PEP_stochastic_perelson.py
# Produces: pep_stochastic_vl_uncertainty.png
#           pep_stochastic_results/*.csv
# No bugs in this script. Run with redirected save paths.
# =============================================================================

def run_fig2():
    t0 = time.time()
    print("\n[2/4] Generating Fig 2 — Stochastic VL uncertainty + CSVs…")

    from PEP_stochastic_perelson import plot_stochastic_analysis, save_results_csv

    save_path = os.path.join(OUTPUT_DIR, 'pep_stochastic_vl_uncertainty.png')
    fig, curves = plot_stochastic_analysis(save_path=save_path)
    plt.close('all')

    save_results_csv(output_dir=STOCHASTIC_CSV_DIR)

    # Copy CSVs to output dir as well for convenience
    import shutil
    for fname in ['pep_stochastic_efficacy_curves.csv',
                  'pep_stochastic_summary_timepoints.csv',
                  'pep_vl_knowledge_premium.csv',
                  'pep_stochastic_ci_width_regimes.csv']:
        src = os.path.join(STOCHASTIC_CSV_DIR, fname)
        dst = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(src):
            shutil.copy2(src, dst)

    print(f"      → {save_path}  ({time.time()-t0:.0f}s)")
    print(f"      → CSVs in {STOCHASTIC_CSV_DIR}/")


# =============================================================================
# STEP 3 — City CSVs: aidsvu_city_profiles.py
# Produces: city_vl_profiles.csv
#           city_pep_efficacy_results.csv
#           city_pep_efficacy_24h_counterfactual.csv
#           city_structural_delay_cost.csv
#
# Bug: standalone main() writes to /mnt/user-data/outputs/ (hardcoded).
# Fix: call the individual functions directly and write to OUTPUT_DIR.
# =============================================================================

def run_city_csvs():
    t0 = time.time()
    print("\n[3/4] Generating city CSVs — AIDSVu 34-city analysis…")

    from aidsvu_city_profiles import (
        build_city_profiles,
        run_city_stratified_analysis,
    )

    # Parse all 34 AIDSVu xlsx files from project root
    city_df = build_city_profiles(data_dir=PROJECT_ROOT)

    p1 = os.path.join(OUTPUT_DIR, 'city_vl_profiles.csv')
    city_df.to_csv(p1, index=False)
    print(f"      → {p1}  ({len(city_df)} cities)")

    # PEP efficacy at city-specific structural delay
    pep_results = run_city_stratified_analysis(
        city_df=city_df,
        pep_delay_hours=None,   # use city-specific delay
        n_simulations=5000
    )
    p2 = os.path.join(OUTPUT_DIR, 'city_pep_efficacy_results.csv')
    pep_results.to_csv(p2, index=False)
    print(f"      → {p2}")

    # Counterfactual at 24h for all cities
    pep_24h = run_city_stratified_analysis(
        city_df=city_df,
        pep_delay_hours=24.0,
        n_simulations=5000
    )
    p3 = os.path.join(OUTPUT_DIR, 'city_pep_efficacy_24h_counterfactual.csv')
    pep_24h.to_csv(p3, index=False)
    print(f"      → {p3}")

    # Structural delay cost (merge actual vs counterfactual)
    merged = pep_results[['city', 'pep_mean_efficacy_pct',
                           'structural_delay_h']].merge(
        pep_24h[['city', 'pep_mean_efficacy_pct']].rename(
            columns={'pep_mean_efficacy_pct': 'pep_24h_counterfactual_pct'}),
        on='city'
    )
    merged['efficacy_lost_to_structural_delay_pp'] = (
        merged['pep_24h_counterfactual_pct'] - merged['pep_mean_efficacy_pct']
    ).round(2)
    merged = merged.sort_values('efficacy_lost_to_structural_delay_pp',
                                ascending=False)
    p4 = os.path.join(OUTPUT_DIR, 'city_structural_delay_cost.csv')
    merged.to_csv(p4, index=False)
    print(f"      → {p4}  ({time.time()-t0:.0f}s)")

    return p2, p3, p4  # return paths for step 4


# =============================================================================
# STEP 4 — Fig 3 + Supp Fig S3: city_stratified_figures.py
#
# BUG 1 (lines 33-36): hardcoded upload paths. Fix by monkey-patching
#   the module's DataFrame globals before the plot code runs.
#
# BUG 2 (lines 105-106): wrong legend labels.
#   Correct thresholds from aidsvu_city_profiles.py:
#     low_barrier      < 8h
#     moderate_barrier 8-18h   ← script says "8-30h"  (WRONG)
#     high_barrier     18-30h  ← script says ">30h"   (WRONG)
#     severe_barrier   >= 30h
#
# Fix: patch the legend_patches list after the module executes its
#   top-level code, before savefig is called. Because city_stratified_figures.py
#   runs plot code at import time (top-level statements), we use
#   importlib + source patching via a temp file.
# =============================================================================

def run_fig3_and_s3(csv_efficacy: str,
                    csv_counterfactual: str,
                    csv_delay_cost: str):
    t0 = time.time()
    print("\n[4/4] Generating Fig 3 + Supp S3 — City-stratified figures…")

    # Read the source script
    src_path = os.path.join(PROJECT_ROOT, 'city_stratified_figures.py')
    with open(src_path, 'r') as f:
        source = f.read()

    # ── Fix 1: replace hardcoded upload paths ──────────────────────────────
    source = source.replace(
        "profiles = pd.read_csv('/mnt/user-data/uploads/1773167806907_city_vl_profiles.csv')",
        f"profiles = pd.read_csv(r'{os.path.join(OUTPUT_DIR, 'city_vl_profiles.csv')}')"
    )
    source = source.replace(
        "efficacy  = pd.read_csv('/mnt/user-data/uploads/1773167806907_city_pep_efficacy_results.csv')",
        f"efficacy  = pd.read_csv(r'{csv_efficacy}')"
    )
    source = source.replace(
        "delay_cost = pd.read_csv('/mnt/user-data/uploads/1773167806907_city_structural_delay_cost.csv')",
        f"delay_cost = pd.read_csv(r'{csv_delay_cost}')"
    )
    source = source.replace(
        "cf24       = pd.read_csv('/mnt/user-data/uploads/1773167806907_city_pep_efficacy_24h_counterfactual.csv')",
        f"cf24       = pd.read_csv(r'{csv_counterfactual}')"
    )

    # ── Fix 2: correct legend labels ──────────────────────────────────────
    source = source.replace(
        "label='Moderate barrier (8\u201330h)'",
        "label='Moderate barrier (8\u201318h)'"
    )
    source = source.replace(
        "label='High barrier (>30h)'",
        "label='High barrier (18\u201330h)'"
    )

    # ── Fix 3: redirect output figure paths ───────────────────────────────
    source = source.replace(
        "fig.savefig('/mnt/user-data/outputs/Fig_CityStratified_PEP.png',",
        f"fig.savefig(r'{os.path.join(OUTPUT_DIR, 'Fig_CityStratified_PEP.png')}', "
    )
    source = source.replace(
        "fig2.savefig('/mnt/user-data/outputs/Fig_CityComparison_Focus.png',",
        f"fig2.savefig(r'{os.path.join(OUTPUT_DIR, 'Fig_CityComparison_Focus.png')}', "
    )

    # ── Write patched source to a temp file and exec it ───────────────────
    tmp_path = os.path.join(OUTPUT_DIR, '_city_stratified_figures_patched.py')
    with open(tmp_path, 'w') as f:
        f.write(source)

    # Execute the patched script in a clean namespace
    exec_globals = {'__file__': tmp_path, '__name__': '__main__'}
    with open(tmp_path, 'r') as f:
        exec(compile(f.read(), tmp_path, 'exec'), exec_globals)

    plt.close('all')
    os.remove(tmp_path)

    fig3_path = os.path.join(OUTPUT_DIR, 'Fig_CityStratified_PEP.png')
    s3_path   = os.path.join(OUTPUT_DIR, 'Fig_CityComparison_Focus.png')
    print(f"      → {fig3_path}")
    print(f"      → {s3_path}  ({time.time()-t0:.0f}s)")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    wall_start = time.time()

    run_fig1()
    run_fig2()
    csv_eff, csv_cf, csv_cost = run_city_csvs()
    run_fig3_and_s3(csv_efficacy=csv_eff,
                    csv_counterfactual=csv_cf,
                    csv_delay_cost=csv_cost)

    print("\n" + "=" * 70)
    print("ALL OUTPUTS COMPLETE")
    print(f"Total time: {time.time()-wall_start:.0f}s")
    print(f"Output dir: {OUTPUT_DIR}")
    print("=" * 70)
    print("\nFigures:")
    for fname in ['pep_parenteral_vl_sweep.png',
                  'pep_stochastic_vl_uncertainty.png',
                  'Fig_CityStratified_PEP.png',
                  'Fig_CityComparison_Focus.png']:
        path = os.path.join(OUTPUT_DIR, fname)
        status = "✓" if os.path.exists(path) else "MISSING"
        print(f"  [{status}] {fname}")

    print("\nCSVs:")
    for fname in ['pep_stochastic_efficacy_curves.csv',
                  'pep_stochastic_summary_timepoints.csv',
                  'pep_vl_knowledge_premium.csv',
                  'pep_stochastic_ci_width_regimes.csv',
                  'city_vl_profiles.csv',
                  'city_pep_efficacy_results.csv',
                  'city_pep_efficacy_24h_counterfactual.csv',
                  'city_structural_delay_cost.csv']:
        path = os.path.join(OUTPUT_DIR, fname)
        status = "✓" if os.path.exists(path) else "MISSING"
        print(f"  [{status}] {fname}")

    print("\nZenodo DOI: https://doi.org/10.5281/zenodo.18746065")
    print("GitHub    : https://github.com/Nyx-Dynamics/Prevention-Theorem")
