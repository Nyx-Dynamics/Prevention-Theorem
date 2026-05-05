"""
Monte Carlo Wrapper for Multiscale Within-Host Model
====================================================
Runs N realizations across multiple V0 values to characterize:

  1. Distribution of T_int (integration completion time) by route
  2. Conditional T_int distribution given non-extinction (relevant population)
  3. t_crit at multiple efficacy thresholds, derived from E_PEP(t)
  4. Route compression ratio as a function of V0 (the manuscript's headline claim)

Output: pandas DataFrame with full simulation results, plus summary statistics.
"""

import numpy as np
import pandas as pd
from multiscale_model import (
    WithinHostParameters, simulate_one_realization
)
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os


def run_mc_at_V0(
    V0: float,
    n_realizations: int = 500,
    base_seed: int = 0,
    params: WithinHostParameters = None,
) -> pd.DataFrame:
    """Run Monte Carlo at a single V0. Returns DataFrame with one row
    per realization."""
    if params is None:
        params = WithinHostParameters()

    rows = []
    # Identify route for barrier delay logic
    route = "mucosal" if V0 < 100 else "parenteral"

    for k in range(n_realizations):
        rng = np.random.default_rng(base_seed * 1_000_000 + k)
        result = simulate_one_realization(params, V0=V0, rng=rng)

        # Apply route-specific barrier delay
        if route == "mucosal":
            T_barrier_h = rng.lognormal(mean=np.log(24), sigma=np.log(1.4))
        else:
            T_barrier_h = 0

        handoff_systemic = (
            result['phase1']['handoff_time'] * 24
            if result['phase1'] and result['phase1']['handoff_time'] else None
        )
        T_int_systemic = result['T_int_hours']

        T_seed_total = (T_barrier_h + handoff_systemic) if handoff_systemic is not None else None
        T_int_total = (T_barrier_h + T_int_systemic) if T_int_systemic is not None else None

        rows.append({
            'realization_idx': k,
            'V0': V0,
            'log10_V0': np.log10(V0) if V0 > 0 else -np.inf,
            'outcome': result['outcome'],
            'extincted': result['outcome'] == 'extinct',
            'integrated': result['outcome'] == 'integrated',
            'T_barrier_h': T_barrier_h,
            'T_seed_systemic': handoff_systemic,
            'T_int_systemic': T_int_systemic,
            'T_seed_total': T_seed_total,
            'T_int_total': T_int_total,
            'T_int_hours': T_int_total, # Maintain backward compatibility
        })
    return pd.DataFrame(rows)


def compute_efficacy_curve(
    df_route: pd.DataFrame,
    eps_max: float = 0.95,
    eps_mid: float = 0.50,
    eps_min: float = 0.0,
    eclipse_floor_hours: float = 21.6,
    time_grid_hours: np.ndarray = None,
) -> pd.DataFrame:
    """
    Compute E_PEP(t) and derived t_crit values from the Monte Carlo results.

    For each realization that integrated, T_int is observed. The seeding time
    T_seed is approximated as the eclipse-phase floor (i.e., the earliest
    point at which any infected cell could begin integration, which is when
    the first eclipse cell becomes productive). This is the natural
    interpretation given the model structure.

    Returns the (mean) efficacy curve and threshold-crossing times at
    eta in {0.50, 0.10, 0.05, 0.01}.
    """
    if time_grid_hours is None:
        time_grid_hours = np.linspace(0, 200, 401)

    # Restrict to realizations that integrated (non-extinct)
    integrated = df_route[df_route['integrated']].copy()
    n_integrated = len(integrated)
    n_total = len(df_route)
    p_extinct = (n_total - n_integrated) / n_total

    if n_integrated == 0:
        return None

    # T_int for each realization (hours)
    T_ints = integrated['T_int_total'].values
    # T_seed approximation: eclipse_floor (earliest possible, conservative)
    # In reality, T_seed precedes T_int by ~tau_eclipse for parenteral and
    # by handoff-time for mucosal (where founder bottleneck adds delay).
    # We'll use the per-realization handoff time as T_seed.
    T_seeds = integrated['T_seed_total'].values

    # P_seed(t) = empirical CDF of T_seed
    P_seed = np.array([np.mean(T_seeds <= t) for t in time_grid_hours])
    # P_int(t) = empirical CDF of T_int
    P_int = np.array([np.mean(T_ints <= t) for t in time_grid_hours])

    # E_PEP(t) for non-extinct realizations
    E_PEP_nonextinct = (
        (1 - P_seed) * eps_max
        + (P_seed - P_int) * eps_mid
        + P_int * eps_min
    )

    # E_PEP(t) marginal over extinction:
    # if extinct, PEP doesn't matter (no infection to prevent), so we
    # report the conditional efficacy for actually-exposed-and-infected
    # cohort. This is the policy-relevant quantity.

    # VL DISTRIBUTION AND KINETIC PARAMETERIZATION:
    #
    # Viral load (VL) distribution parameters are derived from empirical
    # population surveillance data (NHBS/NHAS), with unsuppressed PWID
    # approximated by a log-normal distribution centered at mean log10 ≈ 4.5.
    #
    # Perelson et al. (1996) provides the within-host kinetic parameters
    # used to anchor the temporal structure of the model, including:
    #
    #     - Virion clearance rate (c ≈ 3.07 day⁻¹)
    #     - Eclipse phase duration (~22h)
    #     - Viral generation time (~62h)
    #
    # These kinetic constraints determine the biological timeline of infection
    # progression but do not define the population-level VL distribution.

    # t_crit at various thresholds
    t_crits = {}
    for eta in [0.80, 0.50, 0.10, 0.05, 0.01]:
        below = np.where(E_PEP_nonextinct < eta)[0]
        t_crits[eta] = time_grid_hours[below[0]] if len(below) > 0 else np.nan

    return {
        'time_grid_hours': time_grid_hours,
        'P_seed': P_seed,
        'P_int': P_int,
        'E_PEP': E_PEP_nonextinct,
        't_crits': t_crits,
        'n_integrated': n_integrated,
        'n_total': n_total,
        'p_extinct': p_extinct,
        'T_int_median': np.median(T_ints),
        'T_int_mean': np.mean(T_ints),
        'T_int_p5': np.percentile(T_ints, 5),
        'T_int_p95': np.percentile(T_ints, 95),
    }


def main_route_comparison(
    n_realizations: int = 500,
    V0_grid: list = None,
):
    """Run the route comparison Monte Carlo across multiple V0 values."""
    if V0_grid is None:
        # Spans mucosal (V0=1) through parenteral (V0=10^3) through
        # acute parenteral (V0=10^4)
        V0_grid = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]

    print("="*70)
    print(f"MULTISCALE MONTE CARLO — {n_realizations} realizations per V0")
    print("="*70)

    all_dfs = []
    summaries = []

    for V0 in V0_grid:
        t0 = time.time()
        print(f"\nRunning V0 = {V0:>6.0f} ({n_realizations} realizations)...")
        df = run_mc_at_V0(V0, n_realizations=n_realizations, base_seed=42)
        all_dfs.append(df)
        elapsed = time.time() - t0

        analysis = compute_efficacy_curve(df)
        if analysis is None:
            print(f"  All realizations extinct at V0={V0}")
            continue

        summary = {
            'V0': V0,
            'log10_V0': np.log10(V0),
            'n_total': analysis['n_total'],
            'n_integrated': analysis['n_integrated'],
            'p_extinct': analysis['p_extinct'],
            'T_int_median_h': analysis['T_int_median'],
            'T_int_p5_h': analysis['T_int_p5'],
            'T_int_p95_h': analysis['T_int_p95'],
            't_crit_eta_05_h': analysis['t_crits'][0.05],
            't_crit_eta_50_h': analysis['t_crits'][0.50],
            't_crit_eta_80_h': analysis['t_crits'][0.80],
            'elapsed_sec': elapsed,
        }
        summaries.append(summary)

        print(f"  Extinct: {analysis['p_extinct']*100:>5.1f}%  "
              f"|  T_int median: {analysis['T_int_median']:>5.1f}h  "
              f"(IQR: {analysis['T_int_p5']:.1f}-{analysis['T_int_p95']:.1f})  "
              f"|  t_crit(η=0.05): {analysis['t_crits'][0.05]:>5.1f}h  "
              f"|  Time: {elapsed:.1f}s")

    summary_df = pd.DataFrame(summaries)

    # Combine all per-realization data
    combined = pd.concat(all_dfs, ignore_index=True)

    # Save outputs
    results_dir = 'results/multiscale_model'
    os.makedirs(results_dir, exist_ok=True)
    combined.to_csv(os.path.join(results_dir, 'mc_realizations.csv'), index=False)
    summary_df.to_csv(os.path.join(results_dir, 'mc_summary.csv'), index=False)

    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(summary_df[['V0', 'p_extinct', 'T_int_median_h', 'T_int_p5_h', 'T_int_p95_h',
                      't_crit_eta_05_h', 't_crit_eta_50_h', 't_crit_eta_80_h']].to_string(index=False))

    return summary_df, combined


if __name__ == '__main__':
    # Pilot run with smaller N to confirm everything works,
    # then scale up if the results look right
    summary_df, combined = main_route_comparison(n_realizations=200)
