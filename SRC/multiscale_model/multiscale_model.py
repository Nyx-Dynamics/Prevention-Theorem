"""
Multiscale Within-Host HIV Model
=================================
Replaces the logistic phenomenological model in PEP_parenteral_perelson.py
with a mechanistic multiscale simulation:

  Phase 1 (stochastic): tau-leaping / Gillespie at low founder population
  Phase 2 (deterministic): delay ODE once population is established
  Phase 3 (integration): first-passage time to R(t) >= R*

The eclipse phase is implemented as a discrete delay between cell infection
and integration competence, anchored on Perelson 1996 Table 2 (~22h).
This is the biological floor that gives the integration window its lower
bound regardless of population dynamics.

Author: Built collaboratively for Demidont, AC.
        Designed to replace PEP_parenteral_perelson.py + PEP_stochastic_perelson.py
        with a model whose t_crit values emerge from dynamics rather than
        from hand-set midpoint parameters.
"""

import numpy as np
from scipy.integrate import solve_ivp
from collections import deque
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict
import warnings


# =============================================================================
# PARAMETERS — Perelson 1996 anchored
# =============================================================================
@dataclass
class WithinHostParameters:
    """Mechanistic within-host parameters. All anchored on Perelson 1996
    (Science 271:1582, PMID 8599114) unless explicitly noted as calibrated."""

    # --- Anchored on Perelson 1996 Table 1 ---
    c: float = 23.0          # virion clearance rate (/day)  [Perelson Table 1]
    delta: float = 0.7       # infected cell death rate (/day)  [Perelson Table 1]

    # --- Anchored on Perelson 1996 Table 2 ---
    tau_eclipse: float = 0.9 # eclipse phase (days, ~22h) - MIN time from
                             # cell infection to integration competence
                             # [Perelson Table 2; "minimum duration of
                             # intracellular phase"]

    # --- Standard within-host model values ---
    lam: float = 1e4         # target cell production (cells/mL/day)
    d_T: float = 0.01        # uninfected cell death (/day)
    beta: float = 2.4e-5     # infection rate (mL/virion/day) [standard]
    p: float = 1e3           # virion production (virions/cell/day) [standard]

    # --- Integration kinetics ---
    # alpha is the per-cell-per-day rate that an infected cell, having
    # already passed the eclipse phase, contributes to the latent reservoir.
    # NOTE: this parameter is the one that needs explicit calibration
    # provenance. Either it's anchored on reservoir-seeding kinetics from
    # NHP intravenous challenge (Whitney 2014: reservoir established by
    # day 3) or it's calibrated against another empirical anchor. Whatever
    # value is chosen must be defensible *outside* the CDC 72h window.
    alpha: float = 1e-3      # integration rate post-eclipse (/day)

    # --- Initial conditions ---
    T0: float = 1e6          # initial target cells (CD4 cells/mL peripheral blood)

    # --- Operational thresholds ---
    R_star: float = 10.0     # integration completion threshold (cells/mL)
    I_handoff: float = 10.0  # population at which we hand off from
                             # stochastic to deterministic dynamics
                             # (above this, branching process noise is small)

    def basic_within_host_R0(self) -> float:
        """Within-host R0: expected secondary infected cells per primary
        infected cell. R0 > 1 means infection can establish."""
        # Expected virions produced per infected cell over its lifetime: p/delta
        # Probability each virion infects a target cell: beta*T0 / (c + beta*T0)
        virions_per_cell = self.p / self.delta
        p_infect_per_virion = (self.beta * self.T0) / (self.c + self.beta * self.T0)
        return virions_per_cell * p_infect_per_virion

    def p_single_virion_extinction(self) -> float:
        """Probability that a single founder virion fails to establish
        infection. Combines two failure modes:

        (1) virion clearance before successful target-cell infection,
            P(clear) = c / (c + beta*T0)
        (2) infected cell dies during eclipse phase before becoming productive,
            P(eclipse_death) = 1 - exp(-delta * tau_eclipse)

        These compound multiplicatively: a virion that gets cleared cannot
        infect, and a cell that dies during eclipse cannot establish
        productive infection.
        """
        p_clear = self.c / (self.c + self.beta * self.T0)
        p_eclipse_death = 1.0 - np.exp(-self.delta * self.tau_eclipse)
        # P(single virion fails) = P(clear) + P(infect AND eclipse death)
        return p_clear + (1 - p_clear) * p_eclipse_death


# =============================================================================
# PHASE 1: STOCHASTIC FOUNDER DYNAMICS (tau-leaping)
# =============================================================================
def simulate_founder_phase(
    params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
    dt: float = 0.001,        # 1.4 minute timestep (in days)
    max_time_days: float = 5.0,
) -> Dict:
    """
    Simulate the stochastic founder phase using tau-leaping.

    Tracks individual virions and infected cells until either:
      (a) extinction (V=0 and I=0 with no eclipse cells in pipeline)
      (b) handoff to deterministic phase (I >= I_handoff)
      (c) max_time exceeded

    Returns trajectory and final state for handoff to Phase 2.
    """
    # State: V (free virions), E (eclipse cells, queue of remaining eclipse times),
    #        I (productive infected cells), R (integrated cells)
    V = float(V0)
    E_queue = deque()   # each element is the time-to-end-eclipse for that cell
    I = 0.0
    R = 0.0

    # Trajectory storage (sampled at dt intervals)
    traj_t = [0.0]
    traj_V = [V]
    traj_I = [0.0]
    traj_E = [0.0]
    traj_R = [0.0]

    t = 0.0
    n_steps = int(max_time_days / dt)

    for step in range(n_steps):
        # Compute event rates
        T = params.T0  # target cells effectively constant in founder phase
        rate_V_clear = params.c * V
        rate_V_infect = params.beta * T * V
        rate_I_death = params.delta * I
        rate_I_integrate = params.alpha * I

        # Tau-leaping: number of events in dt
        n_V_clear   = rng.poisson(rate_V_clear * dt) if rate_V_clear > 0 else 0
        n_V_infect  = rng.poisson(rate_V_infect * dt) if rate_V_infect > 0 else 0
        n_I_death   = rng.poisson(rate_I_death * dt) if rate_I_death > 0 else 0
        n_I_integrate = rng.poisson(rate_I_integrate * dt) if rate_I_integrate > 0 else 0
        n_V_produce = rng.poisson(params.p * I * dt) if I > 0 else 0

        # Bound by available populations to avoid negative counts
        n_V_clear   = min(n_V_clear, int(V))
        n_V_infect  = min(n_V_infect, int(V) - n_V_clear)
        n_I_death   = min(n_I_death, int(I))
        n_I_integrate = min(n_I_integrate, int(I) - n_I_death)

        # Update state
        V = V - n_V_clear - n_V_infect + n_V_produce
        # New infections enter eclipse queue with timer = tau_eclipse
        for _ in range(n_V_infect):
            E_queue.append(params.tau_eclipse)
        # Decrement eclipse timers; cells whose eclipse expired become productive
        new_E_queue = deque()
        for time_remaining in E_queue:
            time_remaining -= dt
            if time_remaining <= 0:
                I += 1.0  # cell becomes productive
            else:
                new_E_queue.append(time_remaining)
        E_queue = new_E_queue

        # Cell death and integration come out of productive pool
        I = max(0.0, I - n_I_death - n_I_integrate)
        R = R + n_I_integrate

        t += dt
        traj_t.append(t)
        traj_V.append(V)
        traj_I.append(I)
        traj_E.append(len(E_queue))
        traj_R.append(R)

        # Termination conditions
        if I >= params.I_handoff:
            return {
                'extincted': False,
                'handoff_time': t,
                'handoff_state': {'V': V, 'I': I, 'E': len(E_queue), 'R': R},
                'traj_t': np.array(traj_t),
                'traj_V': np.array(traj_V),
                'traj_I': np.array(traj_I),
                'traj_E': np.array(traj_E),
                'traj_R': np.array(traj_R),
            }
        if V == 0 and I == 0 and len(E_queue) == 0:
            return {
                'extincted': True,
                'handoff_time': None,
                'handoff_state': None,
                'traj_t': np.array(traj_t),
                'traj_V': np.array(traj_V),
                'traj_I': np.array(traj_I),
                'traj_E': np.array(traj_E),
                'traj_R': np.array(traj_R),
            }

    # Reached max_time without handoff or extinction
    return {
        'extincted': False,
        'handoff_time': t,
        'handoff_state': {'V': V, 'I': I, 'E': len(E_queue), 'R': R},
        'traj_t': np.array(traj_t),
        'traj_V': np.array(traj_V),
        'traj_I': np.array(traj_I),
        'traj_E': np.array(traj_E),
        'traj_R': np.array(traj_R),
        'note': 'reached max_time without handoff or extinction'
    }


# =============================================================================
# PHASE 2: DETERMINISTIC ODE WITH ECLIPSE DELAY
# =============================================================================
def deterministic_phase(
    params: WithinHostParameters,
    handoff_state: Dict,
    handoff_time: float,
    R_target: Optional[float] = None,
    max_time_days: float = 7.0,
) -> Dict:
    """
    Integrate the deterministic ODE from the handoff state until either
    R(t) >= R_target or max_time is reached.

    Uses a delay-differential approximation: the eclipse phase is implemented
    by tracking E (eclipse cells) and I (productive cells) as separate
    compartments with first-order transition E -> I at rate 1/tau_eclipse.
    This is a Markovian approximation to the discrete delay.
    """
    if R_target is None:
        R_target = params.R_star

    # State: [T, E, I, V, R]
    T0 = params.T0  # we'll let target cells deplete now that population is large
    y0 = [T0, handoff_state['E'], handoff_state['I'],
          handoff_state['V'], handoff_state['R']]

    eclipse_rate = 1.0 / params.tau_eclipse  # E -> I transition rate

    def odes(t, y):
        T, E, I, V, R = y
        dT = params.lam - params.d_T*T - params.beta*T*V
        dE = params.beta*T*V - eclipse_rate*E
        dI = eclipse_rate*E - params.delta*I
        dV = params.p*I - params.c*V
        dR = params.alpha*I
        return [dT, dE, dI, dV, dR]

    # Event function: stop when R >= R_target
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
        dense_output=True,
        max_step=0.01,
    )

    if sol.t_events[0].size > 0:
        T_int = sol.t_events[0][0]
        reached_target = True
    else:
        T_int = None
        reached_target = False

    return {
        'reached_target': reached_target,
        'T_int_days': T_int,
        'T_int_hours': T_int * 24 if T_int is not None else None,
        'sol': sol,
    }


# =============================================================================
# FULL MULTISCALE SIMULATION (single realization)
# =============================================================================
def simulate_one_realization(
    params: WithinHostParameters,
    V0: float,
    rng: np.random.Generator,
) -> Dict:
    """Run one full multiscale realization and return T_int (hours) or None
    if extinction or max-time reached."""

    phase1 = simulate_founder_phase(params, V0=V0, rng=rng)

    if phase1['extincted']:
        return {
            'T_int_hours': None,
            'outcome': 'extinct',
            'phase1': phase1,
            'phase2': None,
        }

    phase2 = deterministic_phase(
        params,
        handoff_state=phase1['handoff_state'],
        handoff_time=phase1['handoff_time'],
    )

    if phase2['reached_target']:
        return {
            'T_int_hours': phase2['T_int_hours'],  # already in hours
            'outcome': 'integrated',
            'phase1': phase1,
            'phase2': phase2,
        }
    else:
        return {
            'T_int_hours': None,
            'outcome': 'no_integration_in_window',
            'phase1': phase1,
            'phase2': phase2,
        }


# =============================================================================
# QUICK SANITY CHECK
# =============================================================================
if __name__ == '__main__':
    params = WithinHostParameters()

    print("="*70)
    print("MULTISCALE WITHIN-HOST MODEL — SANITY CHECK")
    print("="*70)
    print(f"\nPerelson-anchored within-host R0: {params.basic_within_host_R0():.2f}")
    print(f"  (R0 >> 1 => infection establishes when founder bottleneck cleared)")
    print(f"\nSingle-virion extinction probability: "
          f"{params.p_single_virion_extinction():.3f}")
    print(f"  (mucosal V0=1: ~{params.p_single_virion_extinction()*100:.0f}% "
          f"of exposures fail at founder bottleneck)")

    print(f"\nEclipse phase (Perelson 1996 Table 2): {params.tau_eclipse*24:.1f} h")
    print(f"  This is the biological floor for integration completion time.")
    print(f"  Integration cannot precede ~{params.tau_eclipse*24:.0f}h "
          f"from cell infection regardless of population dynamics.")

    print("\n" + "="*70)
    print("RUNNING THREE PILOT REALIZATIONS")
    print("="*70)

    rng = np.random.default_rng(42)

    # Test parenteral (V0 = 1000)
    print("\nParenteral exposure (V0 = 1000 virions/mL):")
    for i in range(3):
        result = simulate_one_realization(params, V0=1000, rng=rng)
        if result['outcome'] == 'integrated':
            print(f"  Realization {i+1}: T_int = "
                  f"{result['T_int_hours']:.1f} h "
                  f"(handoff at {result['phase1']['handoff_time']*24:.1f}h)")
        else:
            print(f"  Realization {i+1}: {result['outcome']}")

    print("\nMucosal exposure (V0 = 1 virion/mL):")
    for i in range(5):
        result = simulate_one_realization(params, V0=1, rng=rng)
        if result['outcome'] == 'integrated':
            print(f"  Realization {i+1}: T_int = "
                  f"{result['T_int_hours']:.1f} h "
                  f"(handoff at {result['phase1']['handoff_time']*24:.1f}h)")
        else:
            print(f"  Realization {i+1}: {result['outcome']}")
