import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime
from scipy.stats import lognorm


def run_integrated_analysis():
    # --- 1. CONFIGURATION ---
    # Centralized output for your PyCharm project
    TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')
    output_dir = "/Users/acdstudpro/PycharmProjects/Prevention-Theorem/results"
    os.makedirs(output_dir, exist_ok=True)

    n_exposures = 100
    t_initiation = 24.0  # PEP initiation time (hours)
    rng = np.random.RandomState(42)

    print(f"[{TIMESTAMP}] Initializing Integrated Analysis. Outputs: {output_dir}")

    # --- 2. SOURCE VIRAL LOAD & TRANSMISSION RISK ---
    # Log-uniform sampling of source viral load (Vs)
    viral_loads = 10 ** (rng.uniform(3, 6, n_exposures))
    # Transmission probability proportional to log10(Vs)
    p_trans = 0.001 * np.log10(viral_loads)

    # --- 3. MECHANISTIC HITTING TIME & PEP EFFICACY ---
    # Calibrated distributions for seeding (T_seed) and integration (T_int)
    # We apply a 'viral load compression factor' where higher Vs speeds up integration
    vl_compression = 3.0 / np.log10(viral_loads)

    # Hitting times (Mucosal Ref: 20h/48h | Parenteral Ref: 4h/14h)
    t_seed_m = rng.lognormal(mean=np.log(20.0 * vl_compression), sigma=np.log(1.5))
    t_int_m = rng.lognormal(mean=np.log(48.0 * vl_compression), sigma=np.log(1.4))

    # Calculate E_PEP based on the Prevention Theorem
    # E_PEP = (1 - P_seed)*E_max + (P_seed - P_int)*E_mid
    p_seed = np.mean(t_seed_m <= t_initiation)
    p_int = np.mean(t_int_m <= t_initiation)
    e_pep_at_init = (1 - p_seed) * 0.95 + (p_seed - p_int) * 0.50

    # --- 4. POST-EXPOSURE DYNAMICS (ARMA Simulation) ---
    # If PEP fails (Transition to Irreducible State), we simulate seroconversion
    # ramp-up using an ARMA(2,2) model to represent viral load growth
    def simulate_seroconversion(steps=40):
        # Coefficients representing a system with positive growth (unstable/marginal)
        a = [-0.8, -0.4]  # AR parts
        b = [0.2, -0.1]  # MA parts
        y = np.zeros(steps)
        eps = np.random.normal(0, 0.1, steps)

        # Recursive prediction/growth loop
        for t in range(2, steps):
            y[t] = -a[0] * y[t - 1] - a[1] * y[t - 2] + b[0] * eps[t - 1] + b[1] * eps[t - 2] + 0.5  # Add growth bias
        return y

    # --- 5. DATA AGGREGATION & EXPORT ---
    results_df = pd.DataFrame({
        'exposure_id': np.arange(1, n_exposures + 1),
        'source_viral_load': viral_loads,
        'transmission_prob': p_trans,
        'individual_t_int_hours': t_int_m,
        'pep_success_expectation': e_pep_at_init
    })

    csv_path = os.path.join(output_dir, 'integrated_analysis_results.csv')
    results_df.to_csv(csv_path, index=False)
    print(f"[{TIMESTAMP}] Saved: {csv_path}")

    # --- 6. VISUALIZATION ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Plot A: Window Compression vs Viral Load
    ax1.scatter(np.log10(viral_loads), t_int_m, color='#B2182B', alpha=0.6, label='Integration Time ($T_{int}$)')
    ax1.axhline(t_initiation, color='black', linestyle='--', label='PEP Initiation (24h)')
    ax1.set_title("A. Impact of Viral Load on the Finite Prevention Window")
    ax1.set_xlabel("Source Viral Load ($log_{10}$ copies/mL)")
    ax1.set_ylabel("Hours to Irreducible State")
    ax1.legend()
    ax1.grid(alpha=0.2)

    # Plot B: Seroconversion Ramp-up (ARMA Simulation)
    y_growth = simulate_seroconversion()
    ax2.plot(np.arange(len(y_growth)), y_growth, color='#2166AC', linewidth=2)
    ax2.set_title("B. Simulated Post-Integration Viral Ramp-up (ARMA Growth Model)")
    ax2.set_xlabel("Time Steps Post-Failure")
    ax2.set_ylabel("Viral Signal Intensity (Log-Relative)")
    ax2.grid(alpha=0.2)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'integrated_seroconversion_analysis.png')
    plt.savefig(plot_path, dpi=300)
    print(f"[{TIMESTAMP}] Saved: {plot_path}")

    print(f"[{TIMESTAMP}] Analysis Complete. CSV and Plot saved to {output_dir}")


if __name__ == "__main__":
    run_integrated_analysis()