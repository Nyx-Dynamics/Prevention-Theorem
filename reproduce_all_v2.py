#!/usr/bin/env python3
"""
Reproduce all primary figures for "Finite Prevention Windows" (Science Advances v2).
This script executes the model analysis and figure generation in sequence.
"""

import os
import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime

# Project structure
PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / 'SRC'
REVISION_DIR = PROJECT_ROOT / 'v2_revision'

# Create a timestamped run directory
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
RUNS_DIR = REVISION_DIR / 'runs'
FIGURES_DIR = RUNS_DIR / f"run_{timestamp}"

# Ensure figures directory exists
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def run_command(cmd, env_add=None, cwd=None):
    """Run a shell command with specific PYTHONPATH."""
    print(f"\n>>> Running: {' '.join(cmd)}")
    env = os.environ.copy()
    if env_add:
        env['PYTHONPATH'] = env.get('PYTHONPATH', '') + (':' if env.get('PYTHONPATH') else '') + str(env_add)
    
    result = subprocess.run(cmd, env=env, cwd=cwd)
    if result.returncode != 0:
        print(f"Error executing {' '.join(cmd)}")
        sys.exit(result.returncode)

def main():
    print("=" * 80)
    print("V2 REPRODUCTION PIPELINE: PREVENTION THEOREM")
    print("=" * 80)

    # 1. Figure 1: Route-Dependent PEP Efficacy (Multiscale)
    # This uses cached results from SRC/multiscale_model/results_v3/
    # If a fresh run is needed, run SRC/multiscale_model/run_mc_v3.py first.
    print("\n--- Generating Figure 1 (Route Comparison) ---")
    run_command([sys.executable, str(REVISION_DIR / 'make_figure_1_v2.py')])
    
    # Move Figure 1 outputs to the run directory
    # make_figure_1_v2.py usually saves to v2_revision/figures/ by default
    old_fig1_dir = REVISION_DIR / 'figures'
    for ext in ['png', 'pdf']:
        fig1_src = old_fig1_dir / f'Figure_1_Route_Dependent_PEP_Efficacy_Decay.{ext}'
        if fig1_src.exists():
            fig1_dest = FIGURES_DIR / f'Figure_1_Route_Dependent_PEP_Efficacy_Decay.{ext}'
            fig1_src.replace(fig1_dest)
            print(f"Moved Figure 1 to {fig1_dest}")

    # Capture Figure 1 source data (Kinetics results)
    fig1_data_dir = SRC_DIR / 'multiscale_model' / 'results_v3'
    if fig1_data_dir.exists():
        for csv_file in fig1_data_dir.glob('*.csv'):
            csv_dest = FIGURES_DIR / csv_file.name
            shutil.copy2(csv_file, csv_dest)
            print(f"Copied {csv_file.name} to {csv_dest}")

    # 2. Figure 2: Stochastic Efficacy & VL Uncertainty
    print("\n--- Generating Figure 2 (Stochastic VL Uncertainty) ---")
    # This script saves to pep_stochastic_vl_uncertainty.png in current dir by default,
    # or takes an argument. We'll move it after if needed, but let's check its internal save path.
    run_command([sys.executable, str(SRC_DIR / 'perelson/stochastic/PEP_stochastic_perelson.py')], env_add=str(SRC_DIR))
    
    # Move Figure 2 to the figures directory if it was saved in root
    fig2_default = PROJECT_ROOT / 'pep_stochastic_vl_uncertainty.png'
    if fig2_default.exists():
        fig2_target = FIGURES_DIR / 'Figure_2_Stochastic_Efficacy_VL_Uncertainty.png'
        fig2_default.replace(fig2_target)
        print(f"Moved Figure 2 to {fig2_target}")

    # Move Stochastic CSVs to the run directory
    stochastic_csv_dir = PROJECT_ROOT / 'pep_stochastic_results'
    if stochastic_csv_dir.exists():
        for csv_file in stochastic_csv_dir.glob('*.csv'):
            csv_dest = FIGURES_DIR / csv_file.name
            shutil.copy2(csv_file, csv_dest)
            print(f"Copied {csv_file.name} to {csv_dest}")

    # 3. Figure 3 & 4: City-Stratified Analysis
    print("\n--- Generating Figure 3 & 4 (City Stratified) ---")
    run_command([sys.executable, str(REVISION_DIR / 'city_stratified_figures.py')])
    
    # The city_stratified_figures.py script saves to results/city_analysis/
    # We move them to v2_revision/figures/ for the final package
    city_results_dir = PROJECT_ROOT / 'results' / 'city_analysis'

    # Also capture CSVs from city analysis
    for csv_file in city_results_dir.glob('*.csv'):
        csv_dest = FIGURES_DIR / csv_file.name
        shutil.copy2(csv_file, csv_dest)
        print(f"Copied {csv_file.name} to {csv_dest}")

    moves = [
        ('Fig_CityStratified_PEP.png', 'Figure_3_City_Stratified_PEP_Efficacy.png'),
        ('Fig_CityComparison_Focus.png', 'Figure_4_CityComparison_Focus.png'),
        ('Figure_3_City_Stratified_PEP_Efficacy.png', 'Figure_3_City_Stratified_PEP_Efficacy.png'), # handle potential re-runs
        ('Figure_4_CityComparison_Focus.png', 'Figure_4_CityComparison_Focus.png')
    ]
    
    for src_name, dest_name in moves:
        src = city_results_dir / src_name
        if src.exists():
            dest = FIGURES_DIR / dest_name
            src.replace(dest)
            print(f"Moved {src_name} to {dest}")

    print("\n" + "=" * 80)
    print("REPRODUCTION COMPLETE")
    print(f"Figures available in: {FIGURES_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
