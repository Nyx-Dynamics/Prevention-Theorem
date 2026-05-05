"""
Multiscale Within-Host HIV Model (v2 — performance-optimized)
==============================================================
Replaces the queue-based eclipse tracking with population-level dynamics
(O(1) per timestep instead of O(N)). This is biologically equivalent for
large populations (Markovian first-order kinetics) and much faster.

Founder dynamics still tau-leaping; handoff to deterministic ODE happens
when I exceeds threshold; ODE integrates to first-passage of R = R*.

Author: Built collaboratively for Demidont, AC.
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass
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
    """Mechanistic within-host parameters anchored on Perelson 1996."""

    # Anchored on Perelson 1996 Table 1
    c: float = 3.07          # virion clearance rate (/day) [Table 1 mean]
    delta: float = 0.7       # infected cell death rate (/day) [Table 1 mean]

    # Anchored on Perelson 1996 Table 2
    tau_eclipse: float = 0.9 # eclipse phase (days, ~22h) [Table 2]

    # Standard within-host model values
    lam: float = 1e4         # target cell production (cells/mL/day)
    d_T: float = 0.01        # uninfected cell death (/day)
    beta: float = 2.4e-5     # infection rate (mL/virion/day)
    p: float = 1e3           # virion production per cell (virions/cell/day)

    # Integration kinetics
    alpha: float = 1e-3      # integration rate post-eclipse (/day)

    # Initial condition / thresholds
    T0: float = 1e6
    R_star: float = 10.0
    I_handoff: float = 10.0

    def basic_within_host_R0(self) -> float:
        virions_per_cell = self.p / self.delta
        p_infect_per_virion = (self.beta * self.T0) / (self.c + self.beta * self.T0)
        return virions_per_cell * p_infect_per_virion

    def p_single_virion_extinction(self) -> float:
        """Probability single founder virion fails to establish productive
        infection. Two-stage: virion can be cleared, or infected cell can die
        during eclipse before becoming productive."""
        p_clear = self.c / (self.c + self.beta * self.T0)
        p_eclipse_death = 1.0 - np.exp(-self.delta * self.tau_eclipse)
        return p_clear + (1 - p_clear) * p_eclipse_death


# =============================================================================
# PHASE 1: STOCHASTIC FOUNDER DYNAMICS (population-level, tau-leap)
# =============================================================================
def simulate_founder_phase(
    params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
    dt: float = 0.005,        # 7 minute timestep
    max_time_days: float = 5.0,
) -> Dict:
    """
    Tau-leaping with discrete-delay eclipse buffer.

    The eclipse phase is implemented as a fixed-length FIFO buffer: cells
    that get infected at time t become productive exactly at time t+tau_eclipse.
    This is biologically more accurate than a first-order Markov
    approximation (which would give exponentially-distributed eclipse times
    with huge variance).

    Buffer size = tau_eclipse / dt; independent of population size, so
    runtime is O(1) per step regardless of how many cells are in eclipse.
    """
    V = float(V0)
    R = 0.0
    I = 0.0  # productive cells (cells that have exited eclipse)

    # Eclipse buffer: index 0 = cells that just entered eclipse,
    # index n_buffer-1 = cells about to exit eclipse next timestep
    n_buffer = int(round(params.tau_eclipse / dt))
    eclipse_buffer = np.zeros(n_buffer)

    traj_t = [0.0]; traj_V = [V]; traj_E = [0.0]; traj_I = [0.0]; traj_R = [0.0]

    t = 0.0
    n_steps = int(max_time_days / dt)
    record_every = max(1, n_steps // 200)

    for step in range(n_steps):
        T = params.T0  # target cells effectively constant in founder phase

        rate_V_clear   = params.c * V
        rate_V_infect  = params.beta * T * V
        rate_V_produce = params.p * I
        rate_I_death   = params.delta * I
        rate_I_integrate = params.alpha * I
        # Note: cells in eclipse can also die — we model this by applying
        # delta to the entire eclipse pool each step
        rate_E_death = params.delta * eclipse_buffer.sum()

        n_V_clear   = rng.poisson(max(0, rate_V_clear * dt))
        n_V_infect  = rng.poisson(max(0, rate_V_infect * dt))
        n_V_produce = rng.poisson(max(0, rate_V_produce * dt))
        n_I_death   = rng.poisson(max(0, rate_I_death * dt))
        n_I_integrate = rng.poisson(max(0, rate_I_integrate * dt))

        # Bound to physical populations
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

        # Cells exiting eclipse this step (those at end of buffer)
        n_eclipse_exit = eclipse_buffer[-1]

        # Eclipse death (apply uniformly across the eclipse buffer)
        E_total = eclipse_buffer.sum()
        if E_total > 0:
            n_E_death_total = rng.poisson(max(0, params.delta * E_total * dt))
            n_E_death_total = min(n_E_death_total, E_total)
            # Distribute deaths proportionally across the buffer
            death_fractions = eclipse_buffer / E_total
            eclipse_deaths = death_fractions * n_E_death_total
            eclipse_buffer = np.maximum(0, eclipse_buffer - eclipse_deaths)

        # Shift buffer (cells age by dt)
        eclipse_buffer = np.roll(eclipse_buffer, 1)
        eclipse_buffer[0] = n_V_infect  # newly infected enter the buffer
        # n_eclipse_exit was at end of buffer before shift — they become productive

        # Update populations
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


# =============================================================================
# PHASE 2: DETERMINISTIC ODE
# =============================================================================
def deterministic_phase(
    params: WithinHostParameters,
    handoff_state: Dict,
    handoff_time: float,
    R_target: Optional[float] = None,
    max_time_days: float = 7.0,
) -> Dict:
    """Integrate deterministic ODE with E (eclipse) and I (productive)
    compartments until R = R_target or max_time."""
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
# WRAPPER
# =============================================================================
def simulate_one_realization(
    params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
) -> Dict:
    phase1 = simulate_founder_phase(params, V0=V0, rng=rng)
    if phase1['extincted']:
        return {'T_int_hours': None, 'outcome': 'extinct',
                'phase1': phase1, 'phase2': None}
    phase2 = deterministic_phase(params, phase1['handoff_state'],
                                  phase1['handoff_time'])
    if phase2['reached_target']:
        return {'T_int_hours': phase2['T_int_hours'], 'outcome': 'integrated',
                'phase1': phase1, 'phase2': phase2}
    return {'T_int_hours': None, 'outcome': 'no_integration',
            'phase1': phase1, 'phase2': phase2}


# =============================================================================
# QUICK CHECK
# =============================================================================
if __name__ == '__main__':
    import time
    params = WithinHostParameters()

    print("="*70)
    print("MULTISCALE MODEL v2 — TIMING CHECK")
    print("="*70)

    for V0 in [1, 10, 100, 1000, 10000]:
        rng = np.random.default_rng(42)
        t0 = time.time()
        outcomes = []
        T_ints = []
        for k in range(20):
            rng = np.random.default_rng(k)
            r = simulate_one_realization(params, V0=V0, rng=rng)
            outcomes.append(r['outcome'])
            if r['T_int_hours'] is not None:
                T_ints.append(r['T_int_hours'])
        elapsed = time.time() - t0
        n_extinct = sum(1 for o in outcomes if o == 'extinct')
        n_integrated = sum(1 for o in outcomes if o == 'integrated')
        median_t = np.median(T_ints) if T_ints else float('nan')
        print(f"V0={V0:>5.0f}: {elapsed:>5.2f}s for 20 realizations  |  "
              f"extinct: {n_extinct}/20  |  integrated: {n_integrated}/20  |  "
              f"median T_int: {median_t:.1f}h")
