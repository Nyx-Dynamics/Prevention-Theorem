"""
Multiscale Within-Host Model v3
===============================
Extension of v2 with inter-individual parameter heterogeneity.

Each Monte Carlo realization samples its own kinetic parameters from
log-normal distributions centered on the Perelson 1996 anchored values,
representing host-to-host variability observed in NHP challenge studies
and in human within-host kinetics across cohorts.

Parameters with heterogeneity (CV defaults from PEP_stochastic_perelson.py):
    beta   : infection rate (CV = 0.3)
    c      : virion clearance (CV = 0.3)
    delta  : infected cell death (CV = 0.3)
    alpha  : integration rate (CV = 0.3)

Parameters held fixed across realizations:
    tau_eclipse : eclipse phase duration (Perelson 1996 floor; biological
                  constraint, not subject to inter-individual variation in
                  the sense of CV=0.3 — molecular kinetics of integration
                  competence)
    p, lam, d_T : standard within-host literature, low variability
    T0, R_star, I_handoff : initial conditions and operational thresholds

The point of this extension is to widen the T_int distribution to match
NHP cohort variability, so that single-time-point protection statistics
(e.g., Tsai 1998: 50% protection at PEP delay = 48h IV) are reproducible.

Author: Built collaboratively for Demidont, AC.
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, replace
from typing import Optional, Dict


# =============================================================================
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
# =============================================================================
@dataclass
class WithinHostParameters:
    c: float = 3.07          # virion clearance rate (/day) [Table 1 mean]
    delta: float = 0.7       # infected cell death rate (/day) [Table 1 mean]
    tau_eclipse: float = 0.9 # eclipse phase (days, ~22h) [Table 2]
    lam: float = 1e4
    d_T: float = 0.01
    beta: float = 2.4e-5
    p: float = 1e3
    alpha: float = 1e-3
    T0: float = 1e6
    R_star: float = 10.0
    I_handoff: float = 10.0


# =============================================================================
# HETEROGENEITY SAMPLING
# =============================================================================
def sample_heterogeneous_params(
    base: WithinHostParameters,
    rng: np.random.Generator,
    cv: float = 0.3,
    params_to_perturb: tuple = ('beta', 'c', 'delta', 'alpha'),
) -> WithinHostParameters:
    """
    Sample a parameter set for one Monte Carlo realization.

    Each perturbed parameter X is drawn as X = X0 * exp(Normal(0, sigma^2))
    where sigma = sqrt(log(1 + cv^2)) so that the resulting log-normal
    distribution has mean = X0 and coefficient of variation = cv.

    Parameters not in `params_to_perturb` are held at their base values.
    """
    sigma = np.sqrt(np.log(1.0 + cv**2))
    sample = {}
    for name in params_to_perturb:
        x0 = getattr(base, name)
        sample[name] = float(x0 * np.exp(rng.normal(0.0, sigma)))
    return replace(base, **sample)


# =============================================================================
# PHASE 1 + PHASE 2 — copied from v2 (no changes; the heterogeneity is
# applied to the parameters passed in, not to the simulation logic)
# =============================================================================
def simulate_founder_phase(
    params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
    dt: float = 0.005,
    max_time_days: float = 5.0,
) -> Dict:
    V = float(V0)
    R = 0.0
    I = 0.0
    n_buffer = int(round(params.tau_eclipse / dt))
    eclipse_buffer = np.zeros(n_buffer)

    traj_t = [0.0]; traj_V = [V]; traj_E = [0.0]; traj_I = [0.0]; traj_R = [0.0]
    t = 0.0
    n_steps = int(max_time_days / dt)
    record_every = max(1, n_steps // 200)

    for step in range(n_steps):
        T = params.T0
        rate_V_clear   = params.c * V
        rate_V_infect  = params.beta * T * V
        rate_V_produce = params.p * I
        rate_I_death   = params.delta * I
        rate_I_integrate = params.alpha * I

        n_V_clear   = rng.poisson(max(0, rate_V_clear * dt))
        n_V_infect  = rng.poisson(max(0, rate_V_infect * dt))
        n_V_produce = rng.poisson(max(0, rate_V_produce * dt))
        n_I_death   = rng.poisson(max(0, rate_I_death * dt))
        n_I_integrate = rng.poisson(max(0, rate_I_integrate * dt))

        V_lost = n_V_clear + n_V_infect
        if V_lost > V:
            scale = V / V_lost if V_lost > 0 else 0
            n_V_clear = int(n_V_clear * scale)
            n_V_infect = int(n_V_infect * scale)
        I_lost = n_I_death + n_I_integrate
        if I_lost > I:
            scale = I / I_lost if I_lost > 0 else 0
            n_I_death = int(n_I_death * scale)
            n_I_integrate = int(n_I_integrate * scale)

        n_eclipse_exit = eclipse_buffer[-1]

        E_total = eclipse_buffer.sum()
        if E_total > 0:
            n_E_death_total = rng.poisson(max(0, params.delta * E_total * dt))
            n_E_death_total = min(n_E_death_total, E_total)
            death_fractions = eclipse_buffer / E_total
            eclipse_deaths = death_fractions * n_E_death_total
            eclipse_buffer = np.maximum(0, eclipse_buffer - eclipse_deaths)

        eclipse_buffer = np.roll(eclipse_buffer, 1)
        eclipse_buffer[0] = n_V_infect

        V = max(0.0, V - n_V_clear - n_V_infect + n_V_produce)
        I = max(0.0, I + n_eclipse_exit - n_I_death - n_I_integrate)
        R = R + n_I_integrate

        t += dt
        if step % record_every == 0:
            traj_t.append(t); traj_V.append(V)
            traj_E.append(eclipse_buffer.sum())
            traj_I.append(I); traj_R.append(R)

        if I >= params.I_handoff:
            traj_t.append(t); traj_V.append(V)
            traj_E.append(eclipse_buffer.sum()); traj_I.append(I); traj_R.append(R)
            return {
                'extincted': False,
                'handoff_time': t,
                'handoff_state': {'V': V, 'I': I,
                                   'E': eclipse_buffer.sum(),
                                   'eclipse_buffer': eclipse_buffer.copy(),
                                   'R': R},
                'traj_t': np.array(traj_t),
                'traj_V': np.array(traj_V),
                'traj_E': np.array(traj_E),
                'traj_I': np.array(traj_I),
                'traj_R': np.array(traj_R),
            }
        if V == 0 and I == 0 and eclipse_buffer.sum() == 0:
            return {
                'extincted': True,
                'handoff_time': None,
                'handoff_state': None,
                'traj_t': np.array(traj_t),
                'traj_V': np.array(traj_V),
                'traj_E': np.array(traj_E),
                'traj_I': np.array(traj_I),
                'traj_R': np.array(traj_R),
            }

    return {
        'extincted': False,
        'handoff_time': t,
        'handoff_state': {'V': V, 'I': I,
                           'E': eclipse_buffer.sum(),
                           'eclipse_buffer': eclipse_buffer.copy(),
                           'R': R},
        'traj_t': np.array(traj_t),
        'traj_V': np.array(traj_V),
        'traj_E': np.array(traj_E),
        'traj_I': np.array(traj_I),
        'traj_R': np.array(traj_R),
        'note': 'reached max_time'
    }


def deterministic_phase(
    params: WithinHostParameters,
    handoff_state: Dict,
    handoff_time: float,
    R_target: Optional[float] = None,
    max_time_days: float = 7.0,
) -> Dict:
    if R_target is None:
        R_target = params.R_star

    eclipse_rate = 1.0 / params.tau_eclipse
    y0 = [params.T0, handoff_state['E'], handoff_state['I'],
          handoff_state['V'], handoff_state['R']]

    def odes(t, y):
        T, E, I, V, R = y
        dT = params.lam - params.d_T*T - params.beta*T*V
        dE = params.beta*T*V - eclipse_rate*E
        dI = eclipse_rate*E - params.delta*I
        dV = params.p*I - params.c*V
        dR = params.alpha*I
        return [dT, dE, dI, dV, dR]

    def hit_R_target(t, y):
        return y[4] - R_target
    hit_R_target.terminal = True
    hit_R_target.direction = 1

    sol = solve_ivp(
        odes,
        t_span=(handoff_time, handoff_time + max_time_days),
        y0=y0,
        method='LSODA',
        events=hit_R_target,
        rtol=1e-6, atol=1e-9,
        max_step=0.005,
    )

    if sol.t_events[0].size > 0:
        T_int_days = sol.t_events[0][0]
        return {'reached_target': True, 'T_int_days': T_int_days,
                'T_int_hours': T_int_days * 24, 'sol': sol}
    return {'reached_target': False, 'T_int_days': None,
            'T_int_hours': None, 'sol': sol}


# =============================================================================
# WRAPPER WITH HETEROGENEITY
# =============================================================================
def simulate_one_realization_heterogeneous(
    base_params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
    cv: float = 0.3,
    params_to_perturb: tuple = ('beta', 'c', 'delta', 'alpha'),
) -> Dict:
    """Run one realization with sampled (heterogeneous) parameters."""
    sampled = sample_heterogeneous_params(base_params, rng, cv, params_to_perturb)

    phase1 = simulate_founder_phase(sampled, V0=V0, rng=rng)
    if phase1['extincted']:
        return {'T_int_hours': None, 'outcome': 'extinct',
                'phase1': phase1, 'phase2': None,
                'sampled_params': {k: getattr(sampled, k)
                                    for k in params_to_perturb}}

    phase2 = deterministic_phase(sampled, phase1['handoff_state'],
                                  phase1['handoff_time'])
    if phase2['reached_target']:
        return {'T_int_hours': phase2['T_int_hours'],
                'outcome': 'integrated',
                'phase1': phase1, 'phase2': phase2,
                'sampled_params': {k: getattr(sampled, k)
                                    for k in params_to_perturb}}
    return {'T_int_hours': None, 'outcome': 'no_integration',
            'phase1': phase1, 'phase2': phase2,
            'sampled_params': {k: getattr(sampled, k)
                                for k in params_to_perturb}}


# Alias for the noise-free version (CV=0)
def simulate_one_realization_noisefree(
    base_params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
) -> Dict:
    """Convenience: runs with no parameter heterogeneity (matches v2)."""
    return simulate_one_realization_heterogeneous(
        base_params, V0, rng, cv=0.0,
        params_to_perturb=())


# =============================================================================
# QUICK SANITY CHECK
# =============================================================================
if __name__ == '__main__':
    import time
    base = WithinHostParameters()
    rng = np.random.default_rng(42)

    print("="*72)
    print("MULTISCALE v3 — INTER-INDIVIDUAL HETEROGENEITY (CV=0.3)")
    print("="*72)

    print("\nSample 5 perturbed parameter sets:")
    print(f"{'real':<5} {'beta':>14} {'c':>10} {'delta':>10} {'alpha':>14}")
    for i in range(5):
        s = sample_heterogeneous_params(base, rng, cv=0.3)
        print(f"{i+1:<5} {s.beta:>14.3e} {s.c:>10.2f} {s.delta:>10.3f} "
              f"{s.alpha:>14.3e}")

    print("\nBaseline (no heterogeneity):")
    print(f"      {base.beta:>14.3e} {base.c:>10.2f} {base.delta:>10.3f} "
          f"{base.alpha:>14.3e}")

    print("\nPilot timing — 30 realizations per V0 with CV=0.3 noise:")
    for V0 in [1, 100, 1000, 10000]:
        T_ints = []
        outcomes = []
        t0 = time.time()
        for k in range(30):
            r = simulate_one_realization_heterogeneous(
                base, V0=V0, rng=np.random.default_rng(k), cv=0.3
            )
            outcomes.append(r['outcome'])
            if r['T_int_hours']:
                T_ints.append(r['T_int_hours'])
        elapsed = time.time() - t0
        n_ext = sum(1 for o in outcomes if o == 'extinct')
        med = np.median(T_ints) if T_ints else float('nan')
        p5  = np.percentile(T_ints, 5) if T_ints else float('nan')
        p95 = np.percentile(T_ints, 95) if T_ints else float('nan')
        print(f"  V0={V0:>5}: {elapsed:>5.1f}s  ext={n_ext}/30  "
              f"T_int median={med:>5.1f}h  [p5-p95: {p5:.1f}-{p95:.1f}]")
