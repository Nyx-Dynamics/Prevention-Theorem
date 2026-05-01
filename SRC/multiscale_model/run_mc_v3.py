"""
Production Monte Carlo for v3 — sweeps CV (heterogeneity level) to find
the level that best matches the published NHP protection data from
manuscript Table 1.

NHP concordance targets (from Demidont manuscript Table 1):
  Tsai 1998 IV  + 24h delay  + 28-day PEP: 100% protection
  Tsai 1998 IV  + 48h delay  + 28-day PEP:  50% protection
  Otten 2000 IVag + 12h delay: 100% protection
  Otten 2000 IVag + 36h delay: 100% protection
  Otten 2000 IVag + 72h delay:  50% protection

The "protection" rate in NHP studies = fraction of animals that did not
seroconvert. In our model, this maps to E_PEP(delay) — the expected PEP
efficacy at that delay across the heterogeneous cohort.

Strategy: run N=500 realizations at each (V0, CV) combination, compute
E_PEP(t) curve from the resulting T_int distribution, then evaluate at
the NHP-relevant timepoints (24h, 48h, 36h, 72h) and compare to NHP
protection rates.
"""

import numpy as np
import pandas as pd
import os
import time
from multiscale_v3 import (
    WithinHostParameters,
    simulate_one_realization_heterogeneous,
)


def run_mc(
    V0: float,
    cv: float,
    n_realizations: int,
    base_seed: int = 42,
    base_params: WithinHostParameters = None,
):
    """Run N realizations at one (V0, CV) combination."""
    if base_params is None:
        base_params = WithinHostParameters()
    rows = []
    for k in range(n_realizations):
        rng = np.random.default_rng(base_seed * 1_000_000 + k)
        r = simulate_one_realization_heterogeneous(base_params, V0, rng, cv=cv)
        rows.append({
            'realization_idx': k,
            'V0': V0,
            'cv': cv,
            'outcome': r['outcome'],
            'extincted': r['outcome'] == 'extinct',
            'integrated': r['outcome'] == 'integrated',
            'T_int_hours': r['T_int_hours'],
            'handoff_time_hours': (
                r['phase1']['handoff_time'] * 24
                if r['phase1'] and r['phase1']['handoff_time'] else None
            ),
            **{f'sampled_{k}': v for k, v in (r.get('sampled_params') or {}).items()},
        })
    return pd.DataFrame(rows)


def compute_E_PEP_curve(
    df: pd.DataFrame,
    eps_max: float = 0.95,
    eps_mid: float = 0.50,
    eps_min: float = 0.0,
    time_grid: np.ndarray = None,
):
    """Compute E_PEP(t) = (1-P_seed)*eps_max + (P_seed-P_int)*eps_mid + P_int*eps_min."""
    if time_grid is None:
        time_grid = np.linspace(0, 200, 401)

    integrated = df[df['integrated']]
    if len(integrated) == 0:
        return time_grid, np.full_like(time_grid, np.nan)

    T_ints = integrated['T_int_hours'].values
    T_seeds = integrated['handoff_time_hours'].values

    P_seed = np.array([np.mean(T_seeds <= t) for t in time_grid])
    P_int  = np.array([np.mean(T_ints  <= t) for t in time_grid])
    E_PEP = (1-P_seed)*eps_max + (P_seed-P_int)*eps_mid + P_int*eps_min
    return time_grid, E_PEP


def nhp_concordance(df: pd.DataFrame, route_label: str):
    """Compare model E_PEP(t) at NHP timepoints to published protection rates."""
    time_grid, E_PEP = compute_E_PEP_curve(df)

    # NHP timepoints from Tsai 1998 (IV / parenteral) and Otten 2000 (mucosal)
    if 'parenteral' in route_label.lower():
        nhp_data = [
            ('Tsai 1998 IV',   24, 1.00),
            ('Tsai 1998 IV',   48, 0.50),
        ]
    else:
        nhp_data = [
            ('Otten 2000',  12, 1.00),
            ('Otten 2000',  36, 1.00),
            ('Otten 2000',  72, 0.50),
        ]

    rows = []
    for study, delay, nhp_protection in nhp_data:
        idx = np.argmin(np.abs(time_grid - delay))
        model_eff = E_PEP[idx]
        rows.append({
            'study': study,
            'route': route_label,
            'delay_h': delay,
            'NHP_protection': nhp_protection,
            'model_E_PEP': float(model_eff),
            'concordant': abs(model_eff - nhp_protection) < 0.20,  # within 20pp
        })
    return pd.DataFrame(rows)


def t_crit_at(df: pd.DataFrame, eta: float):
    time_grid, E_PEP = compute_E_PEP_curve(df)
    below = np.where(E_PEP < eta)[0]
    return float(time_grid[below[0]]) if len(below) > 0 else float('nan')


def main():
    base = WithinHostParameters()

    # Test 3 CV levels: 0 (no noise, like v2), 0.3 (modest noise),
    # 0.5 (substantial noise — closer to NHP cohort variability)
    cv_grid = [0.0, 0.3, 0.5]
    V0_grid = [1, 1000]    # mucosal vs parenteral PWID
    N = 500                 # 500 realizations per (V0, CV) — tight enough CIs

    print("="*72)
    print(f"MULTISCALE v3 — HETEROGENEITY MONTE CARLO  (N={N} per cell)")
    print("="*72)

    all_results = []
    summary = []
    nhp_rows = []

    for V0 in V0_grid:
        route = 'mucosal' if V0 == 1 else 'parenteral'
        for cv in cv_grid:
            t0 = time.time()
            df = run_mc(V0, cv, n_realizations=N, base_params=base)
            elapsed = time.time() - t0
            all_results.append(df)

            integrated = df[df['integrated']]
            n_int = len(integrated)
            p_extinct = (len(df) - n_int) / len(df)

            if n_int > 0:
                T_med = float(np.median(integrated['T_int_hours']))
                T_p5 = float(np.percentile(integrated['T_int_hours'], 5))
                T_p95 = float(np.percentile(integrated['T_int_hours'], 95))
            else:
                T_med = T_p5 = T_p95 = float('nan')

            t_05 = t_crit_at(df, 0.05)
            t_50 = t_crit_at(df, 0.50)

            summary.append({
                'route': route, 'V0': V0, 'cv': cv,
                'p_extinct': p_extinct, 'n_int': n_int,
                'T_int_median': T_med, 'T_int_p5': T_p5, 'T_int_p95': T_p95,
                't_crit_eta_05': t_05, 't_crit_eta_50': t_50,
                'elapsed': elapsed,
            })

            print(f"  {route:<11} V0={V0:>5} cv={cv:.2f}: "
                  f"ext={p_extinct*100:>4.1f}%  "
                  f"T_int median={T_med:>5.1f}h "
                  f"[{T_p5:.1f}-{T_p95:.1f}]  "
                  f"t_crit_05={t_05:>5.1f}h  ({elapsed:.1f}s)")

            # NHP concordance for this (V0, CV) cell
            nhp = nhp_concordance(df, route)
            nhp['cv'] = cv
            nhp_rows.append(nhp)

    sumdf = pd.DataFrame(summary)
    bigdf = pd.concat(all_results, ignore_index=True)
    nhpdf = pd.concat(nhp_rows, ignore_index=True)

    os.makedirs('results_v3', exist_ok=True)
    sumdf.to_csv('results_v3/heterogeneity_summary.csv', index=False)
    bigdf.to_csv('results_v3/heterogeneity_realizations.csv', index=False)
    nhpdf.to_csv('results_v3/nhp_concordance.csv', index=False)

    print("\n" + "="*72)
    print("NHP CONCORDANCE")
    print("="*72)
    print(nhpdf.to_string(index=False))

    print("\n" + "="*72)
    print("HETEROGENEITY SUMMARY")
    print("="*72)
    cols_show = ['route', 'V0', 'cv', 'p_extinct', 'n_int',
                  'T_int_median', 'T_int_p5', 'T_int_p95',
                  't_crit_eta_05', 't_crit_eta_50']
    print(sumdf[cols_show].to_string(index=False))

    # Headline interpretation
    print("\n" + "="*72)
    print("INTERPRETATION")
    print("="*72)
    print("\nWidening of T_int distribution by heterogeneity level:")
    for V0 in V0_grid:
        sub = sumdf[sumdf['V0']==V0]
        for _, row in sub.iterrows():
            iqr = row['T_int_p95'] - row['T_int_p5']
            print(f"  V0={V0:>5}, cv={row['cv']:.2f}: "
                  f"T_int p5-p95 width = {iqr:.1f}h")

    print("\nNHP concordance (target: |model - NHP| < 0.20):")
    for cv in cv_grid:
        sub = nhpdf[nhpdf['cv']==cv]
        n_concordant = sub['concordant'].sum()
        n_total = len(sub)
        print(f"  cv={cv:.2f}: {n_concordant}/{n_total} timepoints concordant")
        for _, row in sub.iterrows():
            mark = "✓" if row['concordant'] else "✗"
            print(f"    {mark} {row['study']:<14} "
                  f"delay={row['delay_h']:>2}h  "
                  f"NHP={row['NHP_protection']*100:>4.0f}%  "
                  f"model={row['model_E_PEP']*100:>5.1f}%")

    return sumdf, bigdf, nhpdf


if __name__ == '__main__':
    main()
