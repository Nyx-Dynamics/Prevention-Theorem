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

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../core_theorem')))

import numpy as np
import pandas as pd
import time
from absorbing_state import (
    WithinHostParameters,
    simulate_one_realization_heterogeneous,
)
from route_biology import ROUTES


def run_mc(
    V0: float,
    cv: float,
    route: str,
    n_realizations: int,
    base_seed: int = 42,
    base_params: WithinHostParameters = None,
):
    """Run N realizations at one (V0, CV, route) combination."""
    if base_params is None:
        base_params = WithinHostParameters()
    rows = []
    route_params = ROUTES[route]
    for k in range(n_realizations):
        rng = np.random.default_rng(base_seed * 1_000_000 + k)
        r = simulate_one_realization_heterogeneous(base_params, V0, rng, cv=cv)

        # Apply route-specific barrier delay
        if route_params["barrier_delay"]:
            # Logic for mucosal barrier delay (using canonical values from route_biology)
            # Existing code used lognormal(mean=log(24), sigma=log(1.4))
            # We'll use systemic_delay_h as a reference or keep existing distribution for now
            # but acknowledging it's a mucosal route explicitly.
            T_barrier_h = rng.lognormal(mean=np.log(24), sigma=np.log(1.4))
        else:
            T_barrier_h = 0

        handoff_systemic = (
            r['phase1']['handoff_time'] * 24
            if r['phase1'] and r['phase1']['handoff_time'] else None
        )
        T_int_systemic = r['T_int_hours']

        T_seed_total = (T_barrier_h + handoff_systemic) if handoff_systemic is not None else None
        T_int_total = (T_barrier_h + T_int_systemic) if T_int_systemic is not None else None

        rows.append({
            'realization_idx': k,
            'V0': V0,
            'cv': cv,
            'outcome': r['outcome'],
            'extincted': r['outcome'] == 'extinct',
            'integrated': r['outcome'] == 'integrated',
            'T_barrier_h': T_barrier_h,
            'T_seed_systemic': handoff_systemic,
            'T_int_systemic': T_int_systemic,
            'T_seed_total': T_seed_total,
            'T_int_total': T_int_total,
            'T_int_hours': T_int_total, # Maintain backward compatibility
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

    T_ints = integrated['T_int_total'].values
    T_seeds = integrated['T_seed_total'].values

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

    time_grid, E_PEP = compute_E_PEP_curve(df)
    below = np.where(E_PEP < eta)[0]
    return float(time_grid[below[0]]) if len(below) > 0 else float('nan')


def main():
    base = WithinHostParameters()

    # CV grid for heterogeneity levels
    cv_grid = [0.0, 0.3, 0.5]
    
    # Explicit route and V0 definitions
    scenarios = [
        {"route": "mucosal", "V0": 1},
        {"route": "parenteral", "V0": 1000}
    ]
    N = 500                 # 500 realizations per scenario — tight enough CIs

    print("="*72)
    print(f"MULTISCALE v3 — HETEROGENEITY MONTE CARLO  (N={N} per cell)")
    print(f"ROUTE MODEL: mucosal includes anatomical delay = {'yes' if ROUTES['mucosal']['barrier_delay'] else 'no'}")
    print(f"ROUTE MODEL: parenteral bypasses mucosal barrier = {'yes' if not ROUTES['parenteral']['barrier_delay'] else 'no'}")
    print("="*72)

    all_results = []
    summary = []
    nhp_rows = []

    for scenario in scenarios:
        V0 = scenario["V0"]
        route = scenario["route"]
        for cv in cv_grid:
            t0 = time.time()
            df = run_mc(V0, cv, route, n_realizations=N, base_params=base)
            elapsed = time.time() - t0
            all_results.append(df)

            integrated = df[df['integrated']]
            n_int = len(integrated)
            p_extinct = (len(df) - n_int) / len(df)

            if n_int > 0:
                T_med = float(np.median(integrated['T_int_total']))
                T_p5 = float(np.percentile(integrated['T_int_total'], 5))
                T_p95 = float(np.percentile(integrated['T_int_total'], 95))
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
    for scenario in scenarios:
        V0 = scenario["V0"]
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
