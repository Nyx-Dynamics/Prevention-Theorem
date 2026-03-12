import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime


def run_separated_research_flow():
    # --- 1. CONFIGURATION ---
    TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')
    output_dir = "/Users/acdstudpro/PycharmProjects/Prevention-Theorem/results"
    os.makedirs(output_dir, exist_ok=True)

    n_exposures = 100
    t_initiation = 24.0  # PEP initiation time (hours)
    rng = np.random.RandomState(42)

    # Randomly select source viral load (Vs) for 100 exposures
    viral_loads = 10 ** (rng.uniform(3, 6, n_exposures))
    vl_compression = 3.0 / np.log10(viral_loads)

    # Route-Specific Reference Medians
    # Mucosal: T_seed=20h, T_int=48h | Parenteral: T_seed=4h, T_int=14h
    route_configs = {
        'mucosal': {'seed_ref': 20.0, 'int_ref': 48.0, 'color': '#2166AC'},
        'parenteral': {'seed_ref': 4.0, 'int_ref': 14.0, 'color': '#B2182B'}
    }

    print(f"[{TIMESTAMP}] Starting Separated Analysis. Outputs will be saved to: {output_dir}")

    for route_name, params in route_configs.items():
        # --- 2. MECHANISTIC HITTING TIMES ---
        t_seed = rng.lognormal(mean=np.log(params['seed_ref'] * vl_compression), sigma=np.log(1.5), size=n_exposures)
        t_int = rng.lognormal(mean=np.log(params['int_ref'] * vl_compression), sigma=np.log(1.4), size=n_exposures)
        t_seed = np.minimum(t_seed, t_int)

        # --- 3. INDIVIDUAL PEP SUCCESS (Prevention Theorem) ---
        # Success depends on initiation relative to seeding and integration phases
        e_pep = []
        for i in range(n_exposures):
            if t_initiation < t_seed[i]:
                e_pep.append(0.95)  # Preventable State
            elif t_seed[i] <= t_initiation < t_int[i]:
                e_pep.append(0.50)  # Transition Phase
            else:
                e_pep.append(0.05)  # Irreducible State

        # --- 4. EXPORT ROUTE DATA ---
        df = pd.DataFrame({
            'exposure_id': np.arange(1, n_exposures + 1),
            'source_viral_load': viral_loads,
            't_seed_hours': t_seed,
            't_int_hours': t_int,
            'pep_success_expectation': e_pep
        }).sort_values('source_viral_load')

        csv_filename = f'{route_name}_analysis_results.csv'
        csv_path = os.path.join(output_dir, csv_filename)
        df.to_csv(csv_path, index=False)
        print(f"[{TIMESTAMP}] Saved: {csv_path}")

        # --- 5. VISUALIZATION ---
        plt.figure(figsize=(10, 5))
        plt.scatter(np.log10(df['source_viral_load']), df['t_int_hours'], color=params['color'], alpha=0.6,
                    label=f'{route_name.capitalize()} $T_{{int}}$')
        plt.axhline(t_initiation, color='black', linestyle='--', label='PEP Initiation (24h)')
        plt.title(f"Window Compression Analysis: {route_name.capitalize()} Route")
        plt.xlabel("Source Viral Load ($log_{10}$ copies/mL)")
        plt.ylabel("Hours to Irreducibility ($T_{int}$)")
        plt.legend()
        plt.grid(alpha=0.2)

        plot_filename = f'{route_name}_window_analysis.png'
        plot_path = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_path, dpi=300)
        print(f"[{TIMESTAMP}] Saved: {plot_path}")
        plt.close()

        # --- 6. SUMMARY CONSOLE OUTPUT ---
        breakthrough_pct = (df['t_int_hours'] < t_initiation).mean() * 100
        mean_success = df['pep_success_expectation'].mean()
        print(f"\n[{route_name.upper()} RESULTS]")
        print(f"  - Breakthrough % (T_int < 24h): {breakthrough_pct:.1f}%")
        print(f"  - Mean Population Success: {mean_success:.4f}")
        print(f"  - Data saved to: {csv_filename}")


if __name__ == "__main__":
    run_separated_research_flow()