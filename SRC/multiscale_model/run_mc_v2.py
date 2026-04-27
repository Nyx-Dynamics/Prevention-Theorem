"""
Production Monte Carlo for the multiscale within-host model (v2).
Sweeps V0 across 5 logs to characterize:
  - Extinction probability as a function of V0
  - T_int distribution conditional on non-extinction
  - t_crit at multiple efficacy thresholds eta
  - Route compression ratio (mucosal vs parenteral)
"""

import numpy as np
import pandas as pd
import os
import time
from multiscale_v2 import WithinHostParameters, simulate_one_realization


def run_mc_at_V0(
    V0: float,
    n_realizations: int,
    base_seed: int,
    params: WithinHostParameters = None,
):
    if params is None:
        params = WithinHostParameters()
    rows = []
    for k in range(n_realizations):
        rng = np.random.default_rng(base_seed * 1_000_000 + k)
        result = simulate_one_realization(params, V0=V0, rng=rng)
        rows.append({
            'realization_idx': k,
            'V0': V0,
            'log10_V0': np.log10(V0),
            'outcome': result['outcome'],
            'extincted': result['outcome'] == 'extinct',
            'integrated': result['outcome'] == 'integrated',
            'T_int_hours': result['T_int_hours'],
            'handoff_time_hours': (
                result['phase1']['handoff_time'] * 24
                if result['phase1'] and result['phase1']['handoff_time'] else None
            ),
        })
    return pd.DataFrame(rows)


def analyze_route(
    df: pd.DataFrame,
    eps_max: float = 0.95,
    eps_mid: float = 0.50,
    eps_min: float = 0.0,
    time_grid: np.ndarray = None,
):
    """Compute E_PEP(t) and t_crit values from the MC sample."""
    if time_grid is None:
        time_grid = np.linspace(0, 200, 401)

    integrated = df[df['integrated']].copy()
    n_total = len(df)
    n_int = len(integrated)
    p_extinct = (n_total - n_int) / n_total

    if n_int == 0:
        return None

    T_ints = integrated['T_int_hours'].values
    T_seeds = integrated['handoff_time_hours'].values

    P_seed = np.array([np.mean(T_seeds <= t) for t in time_grid])
    P_int  = np.array([np.mean(T_ints  <= t) for t in time_grid])

    E_PEP = (1-P_seed)*eps_max + (P_seed-P_int)*eps_mid + P_int*eps_min

    t_crits = {}
    for eta in [0.80, 0.50, 0.25, 0.10, 0.05, 0.01]:
        below = np.where(E_PEP < eta)[0]
        t_crits[eta] = float(time_grid[below[0]]) if len(below) > 0 else np.nan

    return {
        'time_grid': time_grid,
        'P_seed': P_seed,
        'P_int': P_int,
        'E_PEP': E_PEP,
        't_crits': t_crits,
        'p_extinct': p_extinct,
        'n_total': n_total,
        'n_integrated': n_int,
        'T_int_median': float(np.median(T_ints)),
        'T_int_mean':   float(np.mean(T_ints)),
        'T_int_p5':     float(np.percentile(T_ints, 5)),
        'T_int_p95':    float(np.percentile(T_ints, 95)),
        'T_seed_median': float(np.median(T_seeds)),
    }


def main():
    V0_grid = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
    N = 1000

    print("="*72)
    print(f"MULTISCALE MONTE CARLO — N={N} per V0, V0 grid: {V0_grid}")
    print("="*72)

    summaries = []
    all_dfs = []
    all_curves = {}

    for V0 in V0_grid:
        t0 = time.time()
        df = run_mc_at_V0(V0, n_realizations=N, base_seed=42)
        all_dfs.append(df)
        an = analyze_route(df)
        elapsed = time.time() - t0

        if an is None:
            print(f"V0={V0}: all extincted")
            continue

        all_curves[V0] = an

        s = {
            'V0': V0,
            'log10_V0': np.log10(V0),
            'p_extinct': an['p_extinct'],
            'n_int': an['n_integrated'],
            'T_int_median_h': an['T_int_median'],
            'T_int_p5_h': an['T_int_p5'],
            'T_int_p95_h': an['T_int_p95'],
            'T_seed_median_h': an['T_seed_median'],
            't_crit_eta_80': an['t_crits'][0.80],
            't_crit_eta_50': an['t_crits'][0.50],
            't_crit_eta_10': an['t_crits'][0.10],
            't_crit_eta_05': an['t_crits'][0.05],
            't_crit_eta_01': an['t_crits'][0.01],
            'elapsed': elapsed,
        }
        summaries.append(s)
        print(
            f"V0={V0:>6.0f}: ext={an['p_extinct']*100:>4.1f}%  "
            f"T_int median={an['T_int_median']:>5.1f}h  "
            f"[{an['T_int_p5']:.1f}-{an['T_int_p95']:.1f}]  "
            f"t_crit_05={an['t_crits'][0.05]:>5.1f}h  "
            f"t_crit_50={an['t_crits'][0.50]:>5.1f}h  "
            f"({elapsed:.1f}s)"
        )

    sumdf = pd.DataFrame(summaries)
    bigdf = pd.concat(all_dfs, ignore_index=True)

    os.makedirs('results', exist_ok=True)
    sumdf.to_csv('results/mc_summary.csv', index=False)
    bigdf.to_csv('results/mc_realizations.csv', index=False)

    # Save efficacy curves
    curves_df = pd.DataFrame({'time_h': all_curves[V0_grid[0]]['time_grid']})
    for V0 in V0_grid:
        if V0 in all_curves:
            curves_df[f'E_PEP_V0_{V0}'] = all_curves[V0]['E_PEP']
            curves_df[f'P_int_V0_{V0}'] = all_curves[V0]['P_int']
    curves_df.to_csv('results/mc_efficacy_curves.csv', index=False)

    print("\n" + "="*72)
    print("SUMMARY (saved to results/mc_summary.csv)")
    print("="*72)
    print(sumdf.to_string(index=False))

    # Headline comparison vs manuscript
    print("\n" + "="*72)
    print("HEADLINE: derived t_crit values vs. manuscript claims")
    print("="*72)

    for V0_label, V0 in [('mucosal (V0=1)', 1),
                          ('mucosal-low (V0=10)', 10),
                          ('parenteral (V0=1000)', 1000),
                          ('acute parenteral (V0=10000)', 10000)]:
        if V0 in all_curves:
            an = all_curves[V0]
            print(f"\n  {V0_label}:")
            print(f"    extinction:        {an['p_extinct']*100:.1f}%")
            print(f"    T_int median:      {an['T_int_median']:.1f}h "
                  f"[{an['T_int_p5']:.1f}-{an['T_int_p95']:.1f}]")
            print(f"    t_crit at η=0.80:  {an['t_crits'][0.80]:.1f}h "
                  "(time to drop below 80% efficacy)")
            print(f"    t_crit at η=0.50:  {an['t_crits'][0.50]:.1f}h "
                  "(time to drop below 50%)")
            print(f"    t_crit at η=0.10:  {an['t_crits'][0.10]:.1f}h")
            print(f"    t_crit at η=0.05:  {an['t_crits'][0.05]:.1f}h "
                  "(manuscript-stated threshold)")

    # Compression ratio
    if 1 in all_curves and 1000 in all_curves:
        comp = all_curves[1]['t_crits'][0.05] / all_curves[1000]['t_crits'][0.05]
        print(f"\n  Mucosal/parenteral compression at η=0.05: {comp:.2f}x")
        print(f"  Manuscript-claimed compression: 3.0x")

    return sumdf, all_curves


if __name__ == '__main__':
    sumdf, curves = main()
