#!/usr/bin/env python3
"""
Reproduce All Figures and CSVs for Science Manuscript
=====================================================
Orchestrates the generation of:
1. Fig 1: PEP parenteral VL sweep
2. Fig 2: PEP stochastic VL uncertainty
3. Fig 3 & Supp Fig S3: City-stratified analysis

Usage:
    python reproduce_all.py
"""

import subprocess
import os
import sys
import shutil
from datetime import datetime

TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')

def run_command(cmd, cwd=None):
    print(f"[{datetime.now().strftime('%H:%M')}] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"ERROR: Command failed with return code {result.returncode}")
        sys.exit(1)

def main():
    print("=" * 80)
    print(f"PREVENTION THEOREM REPRODUCTION SUITE | {TIMESTAMP}")
    print("=" * 80)

    SRC = "SRC"
    
    # 1. Script 1 — PEP_parenteral_perelson.py
    run_command([sys.executable, "PEP_parenteral_perelson.py"], cwd=SRC)
    
    # 2. Script 2 — PEP_stochastic_perelson.py
    run_command([sys.executable, "PEP_stochastic_perelson.py"], cwd=SRC)
    
    # 3a. Script 3a — AIDSVu_city_stratified_perelson.py
    run_command([sys.executable, "AIDSVu_city_stratified_perelson.py"], cwd=SRC)
    
    # Copy CSVs from results/city_analysis to SRC/ so city_stratified_figures.py can find them locally
    RESULTS_DIR = os.path.join("results", "city_analysis")
    csv_files = [
        "city_vl_profiles.csv",
        "city_pep_efficacy_results.csv",
        "city_structural_delay_cost.csv",
        "city_pep_efficacy_24h_counterfactual.csv"
    ]
    for f in csv_files:
        src_path = os.path.join(RESULTS_DIR, f)
        if os.path.exists(src_path):
            shutil.copy(src_path, os.path.join(SRC, f))
            print(f"Copied {f} to {SRC}/")

    # 3b. Script 3b — city_stratified_figures.py
    run_command([sys.executable, "city_stratified_figures.py"], cwd=SRC)

    print("\n" + "=" * 80)
    print("REPRODUCTION COMPLETE")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)

if __name__ == "__main__":
    main()
