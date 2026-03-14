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
     Correct per AIDSVu_city_stratified_perelson.py thresholds (<8h, <18h, <30h, >=30h):
       Moderate barrier: 8-18h
       High barrier:     18-30h
  3. Makefile assumes SRC/ subdir — does not exist.
     All scripts run from project root.
  4. AIDSVu_city_stratified_perelson.py writes to /mnt/user-data/outputs/ when run
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

# ── resolve project root ──────────────────────────────────────────────────────
# Script may live in project root OR in a SRC/ subdirectory.
# Source scripts (AIDSVu_city_stratified_perelson.py, city_stratified_figures.py) and
# AIDSVu xlsx files live in the project ROOT, not in SRC/.
# Walk upward until we find AIDSVu_city_stratified_perelson.py.

_script_dir = os.path.dirname(os.path.abspath(__file__))
_candidate  = _script_dir
# Hardcoded project root for robustness if we know it
if _script_dir.endswith('/SRC'):
    PROJECT_ROOT = os.path.dirname(_script_dir)
else:
    PROJECT_ROOT = _script_dir

AIDSVU_DIR = os.path.join(PROJECT_ROOT, 'aidsvu datasets')

sys.path.insert(0, PROJECT_ROOT)           # AIDSVu_city_stratified_perelson, city_stratified…
if _script_dir != PROJECT_ROOT:
    sys.path.insert(0, _script_dir)        # PEP_parenteral_perelson, PEP_stochastic…

# Some scripts might be in SRC
if os.path.exists(os.path.join(PROJECT_ROOT, 'SRC')):
    sys.path.insert(0, os.path.join(PROJECT_ROOT, 'SRC'))

OUTPUT_DIR = os.path.join(_script_dir, 'outputs')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── stochastic CSV subdir (matches PEP_stochastic_perelson.py default) ────────
STOCHASTIC_CSV_DIR = os.path.join(_script_dir, 'pep_stochastic_results')
os.makedirs(STOCHASTIC_CSV_DIR, exist_ok=True)

print("=" * 70)
print("PREVENTION THEOREM — REPRODUCE ALL FIGURES")
print(f"Script dir   : {_script_dir}")
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
# STEP 3 — City CSVs: AIDSVu_city_stratified_perelson.py
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

    from AIDSVu_city_stratified_perelson import (
        build_city_profiles,
        run_city_stratified_analysis,
    )

    # Parse all 34 AIDSVu xlsx files from project root
    city_df = build_city_profiles(data_dir=AIDSVU_DIR)

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
    src_path = os.path.join(PROJECT_ROOT, 'SRC', 'city_stratified_figures.py')
    if not os.path.exists(src_path):
        src_path = os.path.join(PROJECT_ROOT, 'city_stratified_figures.py')
    with open(src_path, 'r') as f:
        source = f.read()

    # ── Fix 1: replace hardcoded upload paths ──────────────────────────────
    source = source.replace(
        "profiles = pd.read_csv(RESULTS_DIR / 'city_vl_profiles.csv')",
        f"profiles = pd.read_csv(r'{os.path.join(OUTPUT_DIR, 'city_vl_profiles.csv')}')"
    )
    source = source.replace(
        "efficacy  = pd.read_csv(RESULTS_DIR / 'city_pep_efficacy_results.csv')",
        f"efficacy  = pd.read_csv(r'{csv_efficacy}')"
    )
    source = source.replace(
        "delay_cost = pd.read_csv(RESULTS_DIR / 'city_structural_delay_cost.csv')",
        f"delay_cost = pd.read_csv(r'{csv_delay_cost}')"
    )
    source = source.replace(
        "cf24       = pd.read_csv(RESULTS_DIR / 'city_pep_efficacy_24h_counterfactual.csv')",
        f"cf24       = pd.read_csv(r'{csv_counterfactual}')"
    )
    source = source.replace(
        "p_focus = profiles[profiles['city'].isin(focus_cities)].copy()",
        f"p_focus = pd.read_csv(r'{os.path.join(OUTPUT_DIR, 'city_vl_profiles.csv')}')\np_focus = p_focus[p_focus['city'].isin(focus_cities)].copy()"
    )

    # ── Fix 2: correct legend labels (if not already corrected in source) ──
    # Note: Source might have \u2013 (en-dash) or - (hyphen)
    source = source.replace(
        "label='Moderate barrier (8\u201330h)'",
        "label='Moderate barrier (8\u201318h)'"
    )
    source = source.replace(
        "label='High barrier (>30h)'",
        "label='High barrier (18\u201330h)'"
    )
    # Also handle hyphens just in case
    source = source.replace(
        "label='Moderate barrier (8-30h)'",
        "label='Moderate barrier (8-18h)'"
    )
    source = source.replace(
        "label='High barrier (>30h)'",
        "label='High barrier (18-30h)'"
    )

    # ── Fix 4: increase padding/margins so city name labels aren't clipped ────
    source = source.replace(
        "plt.tight_layout(pad=4.0)\nfig.subplots_adjust",
        "plt.tight_layout(pad=4.0)\nplt.subplots_adjust(left=0.18, right=0.92, top=0.9, bottom=0.15)\nfig.subplots_adjust"
    )
    source = source.replace(
        "plt.tight_layout(pad=4.0)\nfig2.subplots_adjust",
        "plt.tight_layout(pad=4.0)\nplt.subplots_adjust(top=0.85, bottom=0.15, left=0.1, right=0.95)\nfig2.subplots_adjust"
    )
    source = source.replace(
        "fig.savefig(RESULTS_DIR / 'Fig_CityStratified_PEP.png',",
        f"fig.savefig(r'{os.path.join(OUTPUT_DIR, 'Fig_CityStratified_PEP.png')}', "
    )
    source = source.replace(
        "fig2.savefig(RESULTS_DIR / 'Fig_CityComparison_Focus.png',",
        f"fig2.savefig(r'{os.path.join(OUTPUT_DIR, 'Fig_CityComparison_Focus.png')}', "
    )
    source = source.replace(
        "with open(RESULTS_DIR / 'city_analysis_manuscript_text.txt', 'w') as f:",
        f"with open(r'{os.path.join(OUTPUT_DIR, 'city_analysis_manuscript_text.txt')}', 'w') as f:"
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
# STEP 5 — Panel extraction: crop source figures into submission files
#
# Source figures are 2x2 or 1x3 multi-panel composites.
# Science submission requires each figure and supplement figure as a
# separate file. Panel boundaries were determined by pixel-level whitespace
# scanning of the actual PNGs (see S10.1 in SupplementaryText_v3.docx).
#
# Panel splits are detected automatically by scanning each figure for the
# largest interior whitespace band (the matplotlib subplot gap). This is
# resolution-independent and works regardless of DPI or figure size.
# Empirically verified fallback pixel values for 2000-px-wide renders:
#   par: X=953, Y=810  |  sto: X=955, Y=831
#   cit: X=1022, Y=938 |  foc: X1=654, X2=1328
#
# S10.1 assignment (from SupplementaryText_v3.docx):
#   Fig. 1      → parenteral panels A + D          → Fig_1_submission.png
#   Fig. 2      → stochastic panels A, C, D        → Fig_2_submission.png
#   Fig. 3      → city panels A, B, D              → Fig_3_submission.png
#   Supp. S1    → parenteral panels B + C          → Supp_Fig_S1.png
#   Supp. S2    → stochastic panel B only          → Supp_Fig_S2.png
#   Supp. S3    → city focus all 3 panels          → Supp_Fig_S3.png  (full file)
#   Supp. S4    → city panel C only               → Supp_Fig_S4.png
# =============================================================================

def _crop(img, box):
    """Crop PIL image to box=(left, top, right, bottom) with 2px padding removed."""
    left, top, right, bottom = box
    # Trim 2px from each interior edge to remove any divider bleed
    return img.crop((left, top, right, bottom))


def _hstack(panels, gap_px=20, bg=(255, 255, 255)):
    """Stack PIL images horizontally with a gap."""
    from PIL import Image as PILImage
    total_w = sum(p.width for p in panels) + gap_px * (len(panels) - 1)
    max_h   = max(p.height for p in panels)
    canvas  = PILImage.new('RGB', (total_w, max_h), bg)
    x = 0
    for p in panels:
        canvas.paste(p, (x, 0))
        x += p.width + gap_px
    return canvas


def _vstack(panels, gap_px=20, bg=(255, 255, 255)):
    """Stack PIL images vertically with a gap."""
    from PIL import Image as PILImage
    max_w   = max(p.width for p in panels)
    total_h = sum(p.height for p in panels) + gap_px * (len(panels) - 1)
    canvas  = PILImage.new('RGB', (max_w, total_h), bg)
    y = 0
    for p in panels:
        canvas.paste(p, (0, y))
        y += p.height + gap_px
    return canvas


def run_panel_extraction():
    from PIL import Image as PILImage

    t0 = time.time()
    print("\n[5/5] Extracting panels into submission figures…")

    # ── Source file paths ────────────────────────────────────────────────────
    src_parenteral = os.path.join(OUTPUT_DIR, 'pep_parenteral_vl_sweep.png')
    src_stochastic = os.path.join(OUTPUT_DIR, 'pep_stochastic_vl_uncertainty.png')
    src_city       = os.path.join(OUTPUT_DIR, 'Fig_CityStratified_PEP.png')
    src_focus      = os.path.join(OUTPUT_DIR, 'Fig_CityComparison_Focus.png')

    # Check all source files exist
    missing = [f for f in [src_parenteral, src_stochastic, src_city, src_focus]
               if not os.path.exists(f)]
    if missing:
        print(f"  WARNING: source files not found (steps 1-4 must run first):")
        for f in missing:
            print(f"    MISSING: {f}")
        return

    # ── Load source images ───────────────────────────────────────────────────
    par  = PILImage.open(src_parenteral).convert('RGB')
    sto  = PILImage.open(src_stochastic).convert('RGB')
    cit  = PILImage.open(src_city).convert('RGB')
    foc  = PILImage.open(src_focus).convert('RGB')

    par_w, par_h = par.size
    sto_w, sto_h = sto.size
    cit_w, cit_h = cit.size
    foc_w, foc_h = foc.size

    print(f"  Source figure sizes:")
    print(f"    parenteral : {par_w} x {par_h}")
    print(f"    stochastic : {sto_w} x {sto_h}")
    print(f"    city       : {cit_w} x {cit_h}")
    print(f"    focus      : {foc_w} x {foc_h}")

    # ── Panel split coordinates — whitespace-detected, resolution-adaptive ───
    # Scans each figure for the largest interior whitespace band to find the
    # true panel boundary. Falls back to empirically verified pixel values for
    # 2000px-wide renders if detection fails.
    import numpy as np

    def _find_splits(img_arr, axis, min_band=30, white_thresh=240, white_frac=0.92):
        """
        Return the centre pixel of the largest interior whitespace band along
        `axis` (0 = row bands → Y split, 1 = col bands → X split).
        `min_band` = minimum band width in pixels to count as a real divider.
        """
        # mean across the other spatial axis + all channels
        if axis == 0:
            brightness = np.mean(img_arr > white_thresh, axis=(1, 2))  # per row
        else:
            brightness = np.mean(img_arr > white_thresh, axis=(0, 2))  # per col

        white_mask = brightness > white_frac
        size = len(white_mask)

        # Walk mask to find all contiguous white bands
        bands = []
        in_band = False
        start = 0
        for i in range(size):
            if white_mask[i] and not in_band:
                start = i; in_band = True
            elif not white_mask[i] and in_band:
                bands.append((i - start, start, i - 1))
                in_band = False
        if in_band:
            bands.append((size - start, start, size - 1))

        # Keep interior bands only (ignore leading/trailing whitespace)
        margin = int(size * 0.08)
        interior = [(w, s, e) for w, s, e in bands
                    if s > margin and e < size - margin and w >= min_band]

        if not interior:
            return None
        # Return centre of the widest interior band
        interior.sort(reverse=True)
        _, s, e = interior[0]
        return (s + e) // 2

    def _splits_2x2(img, fallback_x, fallback_y):
        arr = np.array(img)
        x = _find_splits(arr, axis=1)
        y = _find_splits(arr, axis=0)
        x = x if x is not None else fallback_x
        y = y if y is not None else fallback_y
        return x, y

    def _splits_1x3(img, fallback_x1, fallback_x2):
        arr = np.array(img)
        brightness = np.mean(arr > 240, axis=(0, 2))
        white_mask = brightness > 0.92
        size = len(white_mask)
        margin = int(size * 0.05)
        bands = []
        in_band = False
        start = 0
        for i in range(size):
            if white_mask[i] and not in_band:
                start = i; in_band = True
            elif not white_mask[i] and in_band:
                bands.append((i - start, start, i - 1))
                in_band = False
        interior = sorted(
            [(w, s, e) for w, s, e in bands
             if s > margin and e < size - margin and w >= 30],
            reverse=True
        )
        if len(interior) >= 2:
            # Two largest interior column bands → X1 and X2
            centres = sorted([(s + e) // 2 for _, s, e in interior[:2]])
            return centres[0], centres[1]
        return fallback_x1, fallback_x2

    # Empirically verified fallbacks for 2000-px-wide renders
    PAR_X, PAR_Y   = _splits_2x2(par, fallback_x=953,  fallback_y=810)
    STO_X, STO_Y   = _splits_2x2(sto, fallback_x=955,  fallback_y=831)
    CIT_X, CIT_Y   = _splits_2x2(cit, fallback_x=1022, fallback_y=938)
    FOC_X1, FOC_X2 = _splits_1x3(foc, fallback_x1=654, fallback_x2=1328)

    print(f"  Detected split coords:")
    print(f"    PAR  ({par_w}x{par_h}): X={PAR_X}, Y={PAR_Y}")
    print(f"    STO  ({sto_w}x{sto_h}): X={STO_X}, Y={STO_Y}")
    print(f"    CIT  ({cit_w}x{cit_h}): X={CIT_X}, Y={CIT_Y}")
    print(f"    FOC  ({foc_w}x{foc_h}): X1={FOC_X1}, X2={FOC_X2}")

    # ── Crop individual panels ────────────────────────────────────────────
    # Parenteral (A=TL, B=TR, C=BL, D=BR)
    par_A = _crop(par, (0,     0,     PAR_X,   PAR_Y  ))
    par_B = _crop(par, (PAR_X, 0,     par_w,   PAR_Y  ))
    par_C = _crop(par, (0,     PAR_Y, PAR_X,   par_h  ))
    par_D = _crop(par, (PAR_X, PAR_Y, par_w,   par_h  ))

    # Stochastic (A=TL, B=TR, C=BL, D=BR)
    sto_A = _crop(sto, (0,     0,     STO_X,   STO_Y  ))
    sto_B = _crop(sto, (STO_X, 0,     sto_w,   STO_Y  ))
    sto_C = _crop(sto, (0,     STO_Y, STO_X,   sto_h  ))
    sto_D = _crop(sto, (STO_X, STO_Y, sto_w,   sto_h  ))

    # City (A=TL, B=TR, C=BL, D=BR)
    cit_A = _crop(cit, (0,     0,     CIT_X,   CIT_Y  ))
    cit_B = _crop(cit, (CIT_X, 0,     cit_w,   CIT_Y  ))
    cit_C = _crop(cit, (0,     CIT_Y, CIT_X,   cit_h  ))
    cit_D = _crop(cit, (CIT_X, CIT_Y, cit_w,   cit_h  ))

    # City focus (P1, P2, P3 left-to-right)
    foc_1 = _crop(foc, (0,      0, FOC_X1, foc_h))
    foc_2 = _crop(foc, (FOC_X1, 0, FOC_X2, foc_h))
    foc_3 = _crop(foc, (FOC_X2, 0, foc_w,  foc_h))

    bg_white = (255, 255, 255)
    bg_near  = (250, 250, 250)   # city figures use #FAFAFA background

    # ── Assemble and save submission figures ─────────────────────────────
    # Fig. 1 = parenteral A (top row) + D (bottom-right)
    # Layout: A spans full top, D spans bottom-right — reconstruct as 2-panel vertical
    # A = efficacy curves by VL (PRIMARY)
    # D = population efficacy bound / window by route
    fig1 = _vstack([par_A, par_D], gap_px=16, bg=bg_white)
    p = os.path.join(OUTPUT_DIR, 'Fig_1_submission.png')
    fig1.save(p, dpi=(300, 300))
    print(f"      → Fig_1_submission.png  (parenteral panels A + D, {fig1.size})")

    # Fig. 2 = stochastic A + C + D  (Panel B → supplement)
    # A = main efficacy curves with regime shading
    # C = CI width regime plot
    # D = VL knowledge premium
    # Layout: A on top, C and D side by side on bottom
    sto_CD = _hstack([sto_C, sto_D], gap_px=16, bg=bg_white)
    fig2 = _vstack([sto_A, sto_CD], gap_px=16, bg=bg_white)
    p = os.path.join(OUTPUT_DIR, 'Fig_2_submission.png')
    fig2.save(p, dpi=(300, 300))
    print(f"      → Fig_2_submission.png  (stochastic panels A + C + D, {fig2.size})")

    # Fig. 3 = city A + B + D  (Panel C → supplement as S4)
    # A = ranked city efficacy bar chart
    # B = scatter delay vs efficacy
    # D = structural delay components
    cit_AB = _hstack([cit_A, cit_B], gap_px=16, bg=bg_near)
    fig3 = _vstack([cit_AB, cit_D], gap_px=16, bg=bg_near)
    p = os.path.join(OUTPUT_DIR, 'Fig_3_submission.png')
    fig3.save(p, dpi=(300, 300))
    print(f"      → Fig_3_submission.png  (city panels A + B + D, {fig3.size})")

    # Supp. Fig. S1 = parenteral panels B + C (window compression detail)
    s1 = _hstack([par_B, par_C], gap_px=16, bg=bg_white)
    p = os.path.join(OUTPUT_DIR, 'Supp_Fig_S1.png')
    s1.save(p, dpi=(300, 300))
    print(f"      → Supp_Fig_S1.png       (parenteral panels B + C, {s1.size})")

    # Supp. Fig. S2 = stochastic panel B only (VL distribution densities)
    p = os.path.join(OUTPUT_DIR, 'Supp_Fig_S2.png')
    sto_B.save(p, dpi=(300, 300))
    print(f"      → Supp_Fig_S2.png       (stochastic panel B, {sto_B.size})")

    # Supp. Fig. S3 = full city focus figure (all 3 panels — keep full composite)
    import shutil
    p = os.path.join(OUTPUT_DIR, 'Supp_Fig_S3.png')
    shutil.copy2(src_focus, p)
    print(f"      → Supp_Fig_S3.png       (city focus, full composite, {foc.size})")

    # Supp. Fig. S4 = city panel C only (joint VL x structural delay risk surface)
    p = os.path.join(OUTPUT_DIR, 'Supp_Fig_S4.png')
    cit_C.save(p, dpi=(300, 300))
    print(f"      → Supp_Fig_S4.png       (city panel C, {cit_C.size})")

    print(f"\n      Panel extraction complete. ({time.time()-t0:.0f}s)")


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
    run_panel_extraction()

    print("\n" + "=" * 70)
    print("ALL OUTPUTS COMPLETE")
    print(f"Total time: {time.time()-wall_start:.0f}s")
    print(f"Output dir: {OUTPUT_DIR}")
    print("=" * 70)

    print("\nSource figures (4):")
    for fname in ['pep_parenteral_vl_sweep.png',
                  'pep_stochastic_vl_uncertainty.png',
                  'Fig_CityStratified_PEP.png',
                  'Fig_CityComparison_Focus.png']:
        path = os.path.join(OUTPUT_DIR, fname)
        status = "✓" if os.path.exists(path) else "MISSING"
        print(f"  [{status}] {fname}")

    print("\nSubmission figures (7):")
    for fname in ['Fig_1_submission.png',
                  'Fig_2_submission.png',
                  'Fig_3_submission.png',
                  'Supp_Fig_S1.png',
                  'Supp_Fig_S2.png',
                  'Supp_Fig_S3.png',
                  'Supp_Fig_S4.png']:
        path = os.path.join(OUTPUT_DIR, fname)
        status = "✓" if os.path.exists(path) else "MISSING"
        print(f"  [{status}] {fname}")

    print("\nCSVs (8):")
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
