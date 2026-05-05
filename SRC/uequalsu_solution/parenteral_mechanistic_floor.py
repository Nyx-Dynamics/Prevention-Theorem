"""
Mechanistic floor for parenteral P(transmit) — preview of Option C.

Demonstrates how parenteral U=U-equivalent behaviour can emerge from
founder dynamics in the multiscale model, without extrapolating from
sexual-route U=U evidence.

Idea
----
For parenteral inoculation with source plasma VL and inoculum volume V_inoc:

    V0 ~ Poisson(lambda = VL * V_inoc)

Each founder virion fails independently with probability p_fail (from
core_theorem/absorbing_state.py):

    p_fail = p_clear + (1 - p_clear) * p_eclipse_death
    p_clear         = c / (c + beta * T0)
    p_eclipse_death = 1 - exp(-delta * tau_eclipse)

The compound-Poisson extinction probability is:

    P(no establishment | VL) = E[p_fail^V0]
                             = exp(-VL * V_inoc * (1 - p_fail))
                             = exp(-VL * V_inoc * p_success_eff)

So:

    P(transmit | VL, parenteral) = 1 - exp(-VL * V_inoc * p_success_eff)

Calibration gap
---------------
With Perelson defaults, p_success_per_virion ~ 0.47 — an in-vitro-like
value. Empirically, Baggaley 2006 finds per-act needle-share ~0.63% at
"typical" VL, which requires p_success_eff ~1e-5. The 4-order-of-
magnitude gap reflects in-vivo limits the multiscale model omits
(target-cell accessibility, innate restriction, immune clearance,
spatial localization).

We introduce a tissue-access factor kappa, calibrated so that
P(transmit | VL=VL_calibration, V_inoc=V_INOC_PWID) matches Baggaley
0.63%. With kappa fixed, the *shape* of P(transmit) vs VL is then a
mechanistic prediction, and the low-VL behaviour gives the parenteral
analogue of U=U if (and only if) the founder bottleneck dominates.

The whole point of this file: even with Hughes 2012 sexual scaling
swapped out, low VL parenteral exposure has *some* mechanistic basis
for reduced transmission risk. Whether that basis equals "0 below
VL=200" remains open and untestable without empirical needle-share
discordant-couple data.

Cite for methodology:
    Pearson JE, Krapivsky P, Perelson AS. PLoS Comput Biol 2011;
        7:e1001058 (stochastic founder branching theory).
    Conway JM, Coombs D. PLoS Comput Biol 2011;7:e1002033
        (within-host stochastic threshold for HIV establishment).
    Reeves DB et al. PLoS Pathog 2017;13:e1006179
        (Phase 1 stochastic / Phase 2 ODE split, same as ours).
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'core_theorem'))
sys.path.insert(0, HERE)

from absorbing_state import WithinHostParameters
from transmission_probability import p_transmit


# Inoculum volumes by parenteral exposure type (mL).
# Sources:
#   PWID needle: residual blood ~10-30 uL; midpoint 20 uL.
#     Bouvet E. Curr Opin Infect Dis 2009;22:165 (occupational
#     needlestick risk by inoculum).
#   Hollow-bore needlestick: 1 uL midpoint typical of HCW exposure.
#   Solid suture needle: 0.1 uL; very low transfer volume.
#   Blood transfusion: 250 mL standard unit; massive inoculum.
INOCULUM_VOLUME_ML = {
    'pwid_shared_needle':    2.0e-2,
    'needlestick_hollow':    1.0e-3,
    'needlestick_solid':     1.0e-4,
    'blood_transfusion':     2.5e2,
}

# Calibration anchor: Baggaley 2006 per-act needle-share = 0.63% at
# "typical" unsuppressed PWID source VL ~30,000 cp/mL.
VL_CALIBRATION = 30_000
P_TRANSMIT_CALIBRATION = 6.3e-3


def perelson_p_success_per_virion(p: WithinHostParameters = None) -> float:
    """In-vivo-naive per-virion success probability from absorbing_state.

    This is the "ceiling" — assumes target cells are accessible and no
    innate/spatial constraints. Real per-virion success is lower; that
    is what kappa accounts for.
    """
    if p is None:
        p = WithinHostParameters()
    p_clear = p.c / (p.c + p.beta * p.T0)
    p_eclipse_death = 1.0 - np.exp(-p.delta * p.tau_eclipse)
    p_fail = p_clear + (1 - p_clear) * p_eclipse_death
    return 1.0 - p_fail


def calibrate_kappa(exposure_type: str = 'pwid_shared_needle') -> float:
    """Solve for kappa such that the mechanistic P(transmit) matches
    Baggaley 2006 at (VL_CALIBRATION, V_INOC[exposure_type])."""
    v_inoc = INOCULUM_VOLUME_ML[exposure_type]
    p_succ_naive = perelson_p_success_per_virion()

    # 1 - exp(-VL * v_inoc * kappa * p_succ_naive) = P_TRANSMIT_CALIBRATION
    # => kappa = -ln(1 - P) / (VL * v_inoc * p_succ_naive)
    rhs = -np.log(1.0 - P_TRANSMIT_CALIBRATION)
    kappa = rhs / (VL_CALIBRATION * v_inoc * p_succ_naive)
    return float(kappa)


def p_transmit_mechanistic(viral_load,
                           exposure_type: str = 'pwid_shared_needle',
                           kappa: float = None) -> np.ndarray:
    """Mechanistic per-act transmission probability from founder dynamics.

    Returns 1 - exp(-VL * V_inoc * kappa * p_success_per_virion).

    With kappa calibrated to Baggaley 2006, the *shape* of the curve
    in VL is a prediction of the founder-bottleneck framework.
    """
    v_inoc = INOCULUM_VOLUME_ML[exposure_type]
    p_succ_naive = perelson_p_success_per_virion()
    if kappa is None:
        kappa = calibrate_kappa(exposure_type)

    vl = np.asarray(viral_load, dtype=float)
    return 1.0 - np.exp(-vl * v_inoc * kappa * p_succ_naive)


def plot_mechanistic_vs_loglinear():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(
        'Mechanistic founder-bottleneck floor vs Hughes-2012 log-linear\n'
        'Parenteral P(transmit) by source VL — Option C preview',
        fontsize=12, fontweight='bold', y=1.0,
    )

    vl_grid = np.logspace(0, 7, 400)
    kappa = calibrate_kappa('pwid_shared_needle')

    p_mech_pwid = p_transmit_mechanistic(vl_grid, 'pwid_shared_needle')
    p_mech_hollow = p_transmit_mechanistic(vl_grid, 'needlestick_hollow',
                                            kappa=kappa)
    p_loglinear = np.array([p_transmit(v, 'needle_sharing') for v in vl_grid])

    # ----- left: log-log -----
    ax = axes[0]
    ax.loglog(vl_grid, p_mech_pwid, color='crimson', lw=2.5,
              label='Mechanistic (PWID needle, V_inoc=20 uL)')
    ax.loglog(vl_grid, p_mech_hollow, color='darkorange', lw=2,
              label='Mechanistic (HCW needlestick, V_inoc=1 uL)')
    ax.loglog(vl_grid, p_loglinear, color='steelblue', lw=2, ls='--',
              label='Hughes 2012 log-linear (Option B)')

    ax.axvline(200, color='goldenrod', lw=1.5, ls=':', alpha=0.7)
    ax.text(220, 1e-7, 'PARTNER2\nfloor (sexual)',
            fontsize=8, color='goldenrod')
    ax.axhline(P_TRANSMIT_CALIBRATION, color='gray', lw=1, ls=':', alpha=0.5)
    ax.scatter([VL_CALIBRATION], [P_TRANSMIT_CALIBRATION], color='black',
               zorder=5, s=60, label='Baggaley 2006 anchor (0.63%)')

    ax.set_xlabel('Source plasma VL (copies/mL)')
    ax.set_ylabel('P(transmit | exposure)')
    ax.set_title('Log-log: shape difference at low VL', loc='left',
                 fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(1e-8, 1)

    # ----- right: linear at low VL -----
    ax = axes[1]
    vl_low = np.linspace(0, 2000, 400)
    ax.plot(vl_low, p_transmit_mechanistic(vl_low, 'pwid_shared_needle'),
            color='crimson', lw=2.5, label='Mechanistic (PWID)')
    ax.plot(vl_low, [p_transmit(v, 'needle_sharing') for v in vl_low],
            color='steelblue', lw=2, ls='--',
            label='Hughes 2012 log-linear')

    ax.axvline(200, color='goldenrod', lw=1.5, ls=':', alpha=0.7,
               label='PARTNER2 sexual floor')
    ax.set_xlabel('Source plasma VL (copies/mL)')
    ax.set_ylabel('P(transmit | exposure)')
    ax.set_title('Linear zoom: VL < 2,000', loc='left',
                 fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = os.path.join(HERE, 'parenteral_mechanistic_floor.png')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {os.path.relpath(out, PROJECT_ROOT)}")


def report():
    p_succ = perelson_p_success_per_virion()
    kappa = calibrate_kappa('pwid_shared_needle')

    print("=" * 72)
    print("MECHANISTIC FLOOR — DERIVATION FROM FOUNDER DYNAMICS")
    print("=" * 72)
    print(f"\n  p_success_per_virion (Perelson-anchored, naive)  = {p_succ:.3f}")
    print(f"  Tissue-access factor kappa (calibrated to Baggaley) = {kappa:.3e}")
    print(f"  Effective per-virion success in vivo = "
          f"{p_succ * kappa:.3e}")
    print()

    print("P(transmit | VL, exposure) — mechanistic vs log-linear")
    print("-" * 72)
    print(f"  {'VL':>10}   {'mech PWID':>12}   {'mech HCW':>12}   "
          f"{'Hughes Op-B':>12}")
    for vl in [50, 200, 1000, 10_000, 30_000, 100_000, 1_000_000]:
        p_m = p_transmit_mechanistic(vl, 'pwid_shared_needle')
        p_h = p_transmit_mechanistic(vl, 'needlestick_hollow', kappa=kappa)
        p_b = p_transmit(vl, 'needle_sharing')
        print(f"  {vl:>10}   {p_m:>12.3e}   {p_h:>12.3e}   {p_b:>12.3e}")

    print("""
Reading the table
-----------------
The two models give *opposite* tail behaviour at low VL:

  - Hughes 2012 log-linear (Option B) scales as VL^1.531. At low VL it
    drops faster (~34x per log10 VL decade) — an aggressive U=U-like
    decline. This was *measured* in heterosexual partnerships.

  - Mechanistic compound-Poisson founder model scales nearly linearly
    in VL when (VL * V_inoc * kappa * p_success_eff) is small (the
    typical small-x expansion of 1 - exp(-x) ~ x). It drops only ~10x
    per log10 VL decade. This is what you would expect when parenteral
    bypass already supplies a large founder population per unit VL.

At VL=50, mechanistic PWID gives ~1e-5 vs Hughes log-linear ~2e-7 —
the mechanistic model is ~50x *higher*, i.e. *less* U=U-like at low
VL. This is biologically defensible: parenteral inoculation removes
the mucosal-traversal losses that make sexual transmission so VL-
sensitive at low VL.

What this means for the paper
-----------------------------
1. Adopting Option B for parenteral (Hughes-extrapolated) gives the
   *most* U=U-friendly visual; it may overstate parenteral U=U.
2. Adopting this mechanistic model gives a more conservative, less
   U=U-like parenteral curve — closer to the clinical literature's
   acknowledged uncertainty about PWID U=U.
3. Both models converge above VL ~ 1e5 (large founders saturate the
   bottleneck either way), so the policy implications at high VL are
   identical.

This is the long-term Option C path: replace the Option B prefactor
with a mechanistic compound-Poisson founder model + tissue-access
calibration. The current file is the proof-of-concept.
""")


if __name__ == '__main__':
    report()
    plot_mechanistic_vs_loglinear()
