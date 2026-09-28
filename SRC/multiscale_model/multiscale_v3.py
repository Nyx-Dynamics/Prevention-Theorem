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
# PARAMETERS (same as v2)
# =============================================================================
@dataclass
class WithinHostParameters:
    c: float = 23.0
    delta: float = 0.7
    tau_eclipse: float = 0.9
    lam: float = 1e4
    d_T: float = 0.01
    beta: float = 2.4e-5
    p: float = 1e3
    alpha: float = 1e-3
    T0: float = 1e6
    R_star: float = 10.0
    I_handoff: float = 10.0

    # ------------------------------------------------------------------
    # Derived quantities.  I_handoff is NOT a free constant: it is the
    # infected-cell count above which the branching process is safe to
    # replace with its deterministic limit, and that count follows from
    # beta, delta, c and p.  These methods make the derivation checkable
    # against the same parameter set the simulation actually runs on.
    # ------------------------------------------------------------------
    def basic_within_host_R0(self) -> float:
        """Expected secondary infected cells per infected cell.

        Competition form: a virion either infects a target cell (rate
        beta*T0) or is cleared (rate c), so the per-virion infection
        probability is beta*T0 / (c + beta*T0).  Each productive cell
        emits p/delta virions over its lifetime.  This is preferred over
        the textbook beta*T0*p/(delta*c), which assumes c >> beta*T0 —
        an assumption this parameter set does not satisfy (beta*T0 = 24.0
        against c = 23.0).
        """
        p_infect = (self.beta * self.T0) / (self.c + self.beta * self.T0)
        return (self.p / self.delta) * p_infect

    def p_single_virion_extinction(self) -> float:
        """Probability one free virion fails to yield a productive cell.

        Two stages: the virion is cleared before infecting, or the
        infected cell dies during the eclipse phase before becoming
        productive.
        """
        p_clear = self.c / (self.c + self.beta * self.T0)
        p_eclipse_death = 1.0 - np.exp(-self.delta * self.tau_eclipse)
        return p_clear + (1.0 - p_clear) * p_eclipse_death

    def p_extinct_from_I(self, I: float) -> float:
        """P(extinction | I productive infected cells) = (1/R0)^I."""
        return (1.0 / self.basic_within_host_R0()) ** I

    def min_I_handoff(self, tolerance: float = 1e-6) -> int:
        """Smallest infected-cell count with P(extinction) < tolerance."""
        for i in range(1, 1000):
            if self.p_extinct_from_I(i) < tolerance:
                return i
        return 1000

    def handoff_report(self, tolerance: float = 1e-6) -> dict:
        """Everything needed to justify (or reject) the configured value."""
        return {
            'R0_within_host': self.basic_within_host_R0(),
            'q_single_virion_extinction': self.p_single_virion_extinction(),
            'I_handoff_configured': self.I_handoff,
            'I_handoff_required': self.min_I_handoff(tolerance),
            'p_extinct_at_configured': self.p_extinct_from_I(self.I_handoff),
            'tolerance': tolerance,
            'sufficient': self.p_extinct_from_I(self.I_handoff) < tolerance,
            'margin_factor': self.I_handoff / max(self.min_I_handoff(tolerance), 1),
        }


_HANDOFF_CHECK_DONE = False


def validate_handoff_threshold(params: WithinHostParameters,
                               tolerance: float = 1e-6,
                               verbose: bool = True,
                               cv: float = 0.0,
                               quantile: float = 0.01,
                               n_sample: int = 20000,
                               seed: int = 0) -> dict:
    """Assert the configured I_handoff meets the deterministic-limit criterion.

    `params` must be the BASE parameter set, not one heterogeneous draw.
    When cv > 0 the check also evaluates the criterion at an unfavourable
    quantile of the sampled parameter distribution (default: the 1st
    percentile of R0, i.e. the slowest-growing hosts), so the threshold is
    justified for the tail as well as the centre.

    Raises ValueError if the configured threshold leaves a residual
    extinction probability above `tolerance`, so the handoff can never be
    a silently-inherited literal.
    """
    rep = params.handoff_report(tolerance)

    if cv > 0.0:
        rng = np.random.default_rng(seed)
        sigma = np.sqrt(np.log(1.0 + cv ** 2))
        b = params.beta * np.exp(rng.normal(0.0, sigma, n_sample))
        c = params.c * np.exp(rng.normal(0.0, sigma, n_sample))
        d = params.delta * np.exp(rng.normal(0.0, sigma, n_sample))
        R0s = (params.p / d) * (b * params.T0 / (c + b * params.T0))
        R0_q = float(np.quantile(R0s, quantile))
        p_q = (1.0 / R0_q) ** params.I_handoff
        rep.update({'cv': cv, 'quantile': quantile, 'R0_at_quantile': R0_q,
                    'p_extinct_at_quantile': p_q,
                    'sufficient_at_quantile': p_q < tolerance})
        if not rep['sufficient_at_quantile']:
            raise ValueError(
                f"I_handoff = {params.I_handoff:g} is sufficient at the base "
                f"parameters but not at the {quantile:.0%} quantile of the "
                f"cv={cv} distribution: R0 = {R0_q:.2f} gives "
                f"P(extinction) = {p_q:.3e} > {tolerance:.0e}."
            )

    if not rep['sufficient']:
        raise ValueError(
            f"I_handoff = {rep['I_handoff_configured']:g} leaves "
            f"P(extinction) = {rep['p_extinct_at_configured']:.3e}, above the "
            f"tolerance {tolerance:.0e}. Required I >= {rep['I_handoff_required']} "
            f"at R0 = {rep['R0_within_host']:.1f}."
        )
    if verbose:
        msg = (f"[handoff check] base R0 = {rep['R0_within_host']:.1f}; "
               f"q = {rep['q_single_virion_extinction']:.4f}; "
               f"I_handoff = {rep['I_handoff_configured']:g} "
               f"(required >= {rep['I_handoff_required']} at tol "
               f"{tolerance:.0e}; margin {rep['margin_factor']:.1f}x); "
               f"P(ext | I_handoff) = {rep['p_extinct_at_configured']:.2e}")
        if cv > 0.0:
            msg += (f"; at cv={cv} p{quantile:.0%} R0 = {rep['R0_at_quantile']:.1f} "
                    f"-> P(ext) = {rep['p_extinct_at_quantile']:.2e}")
        print(msg)
    return rep


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
    poisson_seed: bool = False,
) -> Dict:
    """Stochastic founder phase (tau-leaping, discrete-delay eclipse buffer).

    V0 is an exact founder count by default.  Set poisson_seed=True to treat
    it as an EXPECTED count and draw n0 ~ Poisson(V0); required when V0 is
    derived empirically (retained volume x source titer), where V0 < 1 means
    "on average fewer than one particle transferred".

    Event counts are rescaled by int() when they exceed V, so a non-integer
    0 < V < 1 fires no events and can never reach an exact V == 0 extinction
    test; such a value previously ran to max_time and returned
    extincted=False.  The extinction condition is therefore V < 1, which is
    identical to V == 0 for integer-valued trajectories.
    """
    if poisson_seed:
        V0 = float(rng.poisson(max(float(V0), 0.0)))
    V = float(V0)
    R = 0.0
    I = 0.0

    if V < 1.0:
        z = np.array([0.0])
        return {
            'extincted': True, 'handoff_time': None, 'handoff_state': None,
            'traj_t': np.array([0.0]), 'traj_V': np.array([V]),
            'traj_E': z, 'traj_I': z, 'traj_R': z,
        }

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
        if V < 1.0 and I == 0 and eclipse_buffer.sum() == 0:
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
    poisson_seed: bool = False,
) -> Dict:
    """Run one realization with sampled (heterogeneous) parameters."""
    global _HANDOFF_CHECK_DONE
    if not _HANDOFF_CHECK_DONE:
        # Validated on the BASE set (plus a tail quantile), never on one draw.
        # Consumes no simulation RNG, so bit-exact reproducibility is unaffected.
        _HANDOFF_CHECK_DONE = True
        validate_handoff_threshold(base_params, cv=cv, verbose=True)

    sampled = sample_heterogeneous_params(base_params, rng, cv, params_to_perturb)

    phase1 = simulate_founder_phase(sampled, V0=V0, rng=rng,
                                    poisson_seed=poisson_seed)
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
    poisson_seed: bool = False,
) -> Dict:
    """Convenience: runs with no parameter heterogeneity (matches v2)."""
    return simulate_one_realization_heterogeneous(
        base_params, V0, rng, cv=0.0,
        params_to_perturb=(), poisson_seed=poisson_seed)


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
