"""
Targeted heterogeneity sweep — try wider noise on alpha and T0 specifically,
since the previous run showed NHP concordance failed not from too-narrow
distributions but from the wrong parameters being perturbed.

Hypothesis: integration time T_int = tau_eclipse + (post-handoff growth time)
The post-handoff time depends mostly on alpha (integration rate) and T0
(target cell pool that determines how fast I builds up). beta, c, delta
are short-timescale parameters that wash out after handoff.

Test combinations:
  Run A: heterogeneity on (alpha) only, CV=1.0 — wide log-normal on alpha
  Run B: heterogeneity on (alpha, T0), each CV=0.5
  Run C: heterogeneity on (alpha, T0, tau_eclipse), wider distributions
"""

import numpy as np
import pandas as pd
import os
import time
from multiscale_v3 import (
    WithinHostParameters,
    simulate_one_realization_heterogeneous,
)


def run_mc_targeted(
    V0: float,
    cv_dict: dict,
    n_realizations: int,
    base_seed: int = 42,
):
    """Run with per-parameter CV values."""
    base = WithinHostParameters()
    rows = []
    for k in range(n_realizations):
        rng = np.random.default_rng(base_seed * 1_000_000 + k)
        # Apply heterogeneity to each parameter individually
        from dataclasses import replace
        sample = {}
        for name, cv in cv_dict.items():
            x0 = getattr(base, name)
            sigma = np.sqrt(np.log(1.0 + cv**2))
            sample[name] = float(x0 * np.exp(rng.normal(0.0, sigma)))
        params_perturbed = replace(base, **sample)

        # Run with no additional heterogeneity (we already perturbed)
        from multiscale_v3 import simulate_founder_phase, deterministic_phase
        phase1 = simulate_founder_phase(params_perturbed, V0=V0, rng=rng)
        if phase1['extincted']:
            outcome = 'extinct'
            T_int_h = None
        else:
            phase2 = deterministic_phase(params_perturbed,
                                          phase1['handoff_state'],
                                          phase1['handoff_time'])
            if phase2['reached_target']:
                outcome = 'integrated'
                T_int_h = phase2['T_int_hours']
            else:
                outcome = 'no_integration'
                T_int_h = None

        rows.append({
            'realization_idx': k,
            'V0': V0,
            'outcome': outcome,
            'extincted': outcome == 'extinct',
            'integrated': outcome == 'integrated',
            'T_int_hours': T_int_h,
            'handoff_time_hours': (
                phase1['handoff_time'] * 24
                if phase1.get('handoff_time') else None
            ),
            **{f'p_{k}': v for k, v in sample.items()},
        })
    return pd.DataFrame(rows)


def E_PEP_at(df, t_query, eps_max=0.95, eps_mid=0.50, eps_min=0.0):
    """Evaluate E_PEP(t) at a specific time."""
    integrated = df[df['integrated']]
    if len(integrated) == 0:
        return float('nan')
    T_ints = integrated['T_int_hours'].values
    T_seeds = integrated['handoff_time_hours'].values
    P_seed = float(np.mean(T_seeds <= t_query))
    P_int  = float(np.mean(T_ints  <= t_query))
    return (1-P_seed)*eps_max + (P_seed-P_int)*eps_mid + P_int*eps_min


def t_crit_at(df, eta, t_max=200):
    times = np.linspace(0, t_max, 401)
    eff = np.array([E_PEP_at(df, t) for t in times])
    below = np.where(eff < eta)[0]
    return float(times[below[0]]) if len(below) > 0 else float('nan')


def assess(df, route_label):
    """Compute summary stats and NHP concordance."""
    integrated = df[df['integrated']]
    n_int = len(integrated)
    p_extinct = (len(df) - n_int) / len(df)

    if n_int > 0:
        T_med = float(np.median(integrated['T_int_hours']))
        T_p5 = float(np.percentile(integrated['T_int_hours'], 5))
        T_p95 = float(np.percentile(integrated['T_int_hours'], 95))
    else:
        T_med = T_p5 = T_p95 = float('nan')

    if 'parenteral' in route_label.lower():
        nhp_pts = [(24, 1.00, 'Tsai 24h'), (48, 0.50, 'Tsai 48h')]
    else:
        nhp_pts = [(12, 1.00, 'Otten 12h'), (36, 1.00, 'Otten 36h'),
                   (72, 0.50, 'Otten 72h')]

    nhp_results = []
    for delay, target, label in nhp_pts:
        model_eff = E_PEP_at(df, delay)
        nhp_results.append({
            'point': label,
            'delay': delay,
            'NHP_target': target,
            'model_E_PEP': model_eff,
            'concordant': abs(model_eff - target) < 0.20,
        })

    return {
        'p_extinct': p_extinct,
        'n_int': n_int,
        'T_int_median': T_med,
        'T_int_p5': T_p5,
        'T_int_p95': T_p95,
        't_crit_05': t_crit_at(df, 0.05),
        't_crit_50': t_crit_at(df, 0.50),
        'nhp': nhp_results,
    }


def run_scenario(label: str, cv_dict: dict, N: int = 500):
    print(f"\n{'='*72}")
    print(f"SCENARIO: {label}")
    print(f"  Heterogeneity: {cv_dict}")
    print('='*72)

    for V0 in [1, 1000]:
        route = 'mucosal' if V0 == 1 else 'parenteral'
        t0 = time.time()
        df = run_mc_targeted(V0, cv_dict, N)
        elapsed = time.time() - t0
        a = assess(df, route)

        print(f"\n  {route:<11} V0={V0:>5} (N={N}, {elapsed:.1f}s)")
        print(f"    extinct: {a['p_extinct']*100:>5.1f}%  "
              f"T_int median: {a['T_int_median']:>5.1f}h  "
              f"[p5-p95: {a['T_int_p5']:.1f}-{a['T_int_p95']:.1f}, "
              f"width {a['T_int_p95']-a['T_int_p5']:.1f}h]")
        print(f"    t_crit at η=0.05: {a['t_crit_05']:>5.1f}h, "
              f"η=0.50: {a['t_crit_50']:>5.1f}h")
        print(f"    NHP concordance:")
        n_concordant = sum(1 for n in a['nhp'] if n['concordant'])
        for n in a['nhp']:
            mark = '✓' if n['concordant'] else '✗'
            print(f"      {mark} {n['point']:<14} "
                  f"NHP={n['NHP_target']*100:>4.0f}%  "
                  f"model={n['model_E_PEP']*100:>5.1f}%")
        print(f"    -> {n_concordant}/{len(a['nhp'])} concordant")


if __name__ == '__main__':
    # Baseline: v3 default (cv=0.3 on β, c, δ, α)
    run_scenario(
        "Baseline v3: CV=0.3 on β, c, δ, α",
        {'beta': 0.3, 'c': 0.3, 'delta': 0.3, 'alpha': 0.3}
    )

    # Run A: alpha-dominated heterogeneity
    run_scenario(
        "A: wide alpha (CV=1.0), modest others",
        {'beta': 0.3, 'c': 0.3, 'delta': 0.3, 'alpha': 1.0}
    )

    # Run B: alpha + T0 heterogeneity
    run_scenario(
        "B: alpha CV=0.7, T0 CV=0.5, others modest",
        {'beta': 0.3, 'c': 0.3, 'delta': 0.3, 'alpha': 0.7, 'T0': 0.5}
    )

    # Run C: alpha + T0 + tau_eclipse heterogeneity
    run_scenario(
        "C: alpha CV=0.7, T0 CV=0.5, tau_eclipse CV=0.3",
        {'beta': 0.3, 'c': 0.3, 'delta': 0.3, 'alpha': 0.7,
         'T0': 0.5, 'tau_eclipse': 0.3}
    )

    print("\n" + "="*72)
    print("DONE — see results above for which heterogeneity profile best")
    print("matches NHP protection data")
    print("="*72)
