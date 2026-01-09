# This script generates boundary-condition data for the Prevention Theorem.
# Post-integration dynamics and reservoir modeling are intentionally excluded
# to preserve the scope of the Epidemics submission (Paper 1).

import numpy as np
import pandas as pd
import os
from PEP_mucosal import InfectionEstablishmentModel

def generate_figure_1_data(output_dir):
    """Generates data for Figure 1: Theoretical Framework from prevention_theorem_figures.py."""
    t = np.linspace(0, 96, 500)  # Hours

    # Sigmoidal functions for biological probabilities
    # Seed: Reservoir seeding starts early
    p_seed = 1 / (1 + np.exp(-0.15 * (t - 24)))
    # Int: Integration lags seeding
    p_int = 1 / (1 + np.exp(-0.2 * (t - 48)))

    # Efficacy Function E_PEP(t)
    # E_max = 0.95, E_mid = 0.5, E_min = 0
    e_max, e_mid, e_min = 0.95, 0.5, 0.0
    e_pep = (1 - p_seed) * e_max + (p_seed - p_int) * e_mid + p_int * e_min

    # Irreversibility Metric (1 - E_PEP)
    irreversibility = 1 - e_pep

    df = pd.DataFrame({
        'hours': t,
        'p_seed': p_seed,
        'p_int': p_int,
        'e_pep': e_pep,
        'irreversibility': irreversibility
    })

    df.to_csv(os.path.join(output_dir, 'fig1_prevention_theorem_data.csv'), index=False)
    print(f"Saved Figure 1 data to {output_dir}")
    return df

def generate_figure_2_data(output_dir):
    """Generates data for Figure 2: Mucosal vs Parenteral from prevention_theorem_figures.py."""
    t = np.linspace(0, 96, 500)

    # Efficacy decay
    # Mucosal: Slow decay (immune bottlenecks)
    eff_mucosal = 1 / (1 + np.exp(0.15 * (t - 60)))
    # Parenteral: Fast decay (systemic access)
    eff_parenteral = 1 / (1 + np.exp(0.3 * (t - 18)))

    df = pd.DataFrame({
        'hours': t,
        'eff_mucosal': eff_mucosal,
        'eff_parenteral': eff_parenteral
    })

    df.to_csv(os.path.join(output_dir, 'fig2_window_compression_data.csv'), index=False)
    print(f"Saved Figure 2 data to {output_dir}")
    return df

def main():
    output_dir = "data_output"
    os.makedirs(output_dir, exist_ok=True)
    
    generate_figure_1_data(output_dir)
    generate_figure_2_data(output_dir)
    
    print(f"\nAll data generated successfully in {output_dir}/")

if __name__ == "__main__":
    main()
