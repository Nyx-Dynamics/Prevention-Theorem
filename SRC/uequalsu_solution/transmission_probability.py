"""
p_transmit(VL, route): per-act transmission probability multiplier
for the Prevention Theorem PEP framework.

This is the missing prefactor that lets the parenteral_route and
mucosal_route phenomenological models respect U=U for sexual exposure
while remaining honest about its absence for parenteral exposure.

Functional form
---------------
    p_transmit(VL, route) = baseline[route] * (VL / VL_anchor)^beta
                            * floor_indicator(VL, route)

The model is the joint prevention quantity:
    P(prevent infection) = p_transmit(VL, route) * P(PEP succeeds | transmit)

At VL<200 sexual, p_transmit = 0 by PARTNER2 — making the joint quantity
exactly zero, which is the correct U=U statement (no transmission to
prevent). For parenteral, no floor is applied: the absence of evidence
is itself the publishable scientific finding.

Sources
-------
Hughes JP et al. JID 2012;205:358 (HPTN 052 secondary analysis).
    Hazard ratio per log10 VL = 2.89 in heterosexual partnerships.
    => beta = log2(2.89) ~ 1.531 used here for both sexual and
       parenteral routes (parenteral choice is an extrapolation,
       documented as such — see PARENTERAL_BETA below).
Patel P et al. AIDS 2014;28:1509. Per-act transmission probability
    estimates by exposure route (Table 1). Anchor for baseline[route].
Rodger AJ et al. Lancet 2019;393:2428 (PARTNER2). Zero phylogenetically
    linked transmissions in ~75,000 condomless acts at VL<200.
    => SEXUAL_VL_FLOOR = 200 cp/mL.
Quinn TC et al. NEJM 2000;342:921 (Rakai). Original log-linear fit;
    no transmissions observed below ~1500 cp/mL across 415 couples.
Cohen MS et al. NEJM 2011;365:493 (HPTN 052 main). 96% reduction in
    linked transmission with early ART — established U=U for sexual.
Bavinton BR et al. Lancet HIV 2018;5:e438 (Opposites Attract). Same
    finding for MSM — confirms PARTNER2 floor extends to anal route.
Baggaley RF et al. AIDS 2006;20:805. Meta-analysis: per-act needle-
    sharing risk = 0.63%. There is no needle-sharing analogue to
    PARTNER. This justifies "U=U has not been validated for
    parenteral exposure" and motivates the asymmetric floor policy.
"""

from typing import Union

import numpy as np


PER_ACT_BASELINE = {
    'receptive_anal':           1.38e-2,
    'insertive_anal':           1.10e-3,
    'receptive_vaginal':        8.00e-4,
    'insertive_vaginal':        4.00e-4,
    'mother_to_child_in_utero': 1.00e-2,
    'needle_sharing':           6.30e-3,
    'needlestick_hollow':       2.30e-3,
    'pwid_shared_needle':       6.30e-3,
}

# Hughes 2012 HPTN 052: HR=2.89 per log10 VL.
# beta = log2(2.89) is the exponent on (VL/VL_anchor).
SEXUAL_BETA = np.log2(2.89)

# Parenteral beta is an extrapolation — Hughes 2012 modeled heterosexual
# transmission only. Two defensible choices:
#   (a) same exponent (mechanistic plausibility, same virus/same biology)
#   (b) weaker exponent (parenteral bypass less VL-sensitive because
#       founder population is larger regardless of source VL)
# Default to (a) but document the assumption; reviewer can challenge.
PARENTERAL_BETA = SEXUAL_BETA

VL_ANCHOR = 10 ** 4.6  # ~40,000 cp/mL — HPTN 052 cohort median

# PARTNER2 floor for sexual routes. Below this VL, observed P(transmit)=0.
SEXUAL_VL_FLOOR = 200  # cp/mL

SEXUAL_ROUTES = {
    'receptive_anal', 'insertive_anal',
    'receptive_vaginal', 'insertive_vaginal',
    'mother_to_child_in_utero',
}
PARENTERAL_ROUTES = {
    'needle_sharing', 'needlestick_hollow', 'pwid_shared_needle',
}


def p_transmit(viral_load: Union[float, np.ndarray],
               route: str,
               apply_uequalsu_floor: bool = True
               ) -> Union[float, np.ndarray]:
    """Per-act transmission probability for a given source VL and route.

    Args:
        viral_load: source plasma VL in copies/mL (scalar or array).
        route: one of PER_ACT_BASELINE keys.
        apply_uequalsu_floor: if True (default), enforce P=0 below
            SEXUAL_VL_FLOOR for sexual routes (PARTNER2). Set False to
            generate the "what if no floor existed" counterfactual.

    Returns:
        Per-act transmission probability in [0, 1]. Scalar in -> scalar
        out; array in -> array out.
    """
    if route not in PER_ACT_BASELINE:
        raise ValueError(
            f"Unknown route '{route}'. Known: {sorted(PER_ACT_BASELINE)}"
        )

    scalar_input = np.ndim(viral_load) == 0
    vl = np.maximum(np.asarray(viral_load, dtype=float), 1.0)

    beta = PARENTERAL_BETA if route in PARENTERAL_ROUTES else SEXUAL_BETA
    p = PER_ACT_BASELINE[route] * (vl / VL_ANCHOR) ** beta

    if apply_uequalsu_floor and route in SEXUAL_ROUTES:
        p = np.where(vl < SEXUAL_VL_FLOOR, 0.0, p)

    p = np.clip(p, 0.0, 1.0)
    return float(p) if scalar_input else p


def joint_prevention(pep_efficacy_given_transmission: Union[float, np.ndarray],
                     viral_load: Union[float, np.ndarray],
                     route: str
                     ) -> Union[float, np.ndarray]:
    """Joint Panel A quantity: P(prevent infection | exposure).

        joint = p_transmit(VL, route) * P(PEP succeeds | transmit)

    At VL<200 sexual, joint = 0 (correct U=U statement). At VL<200
    parenteral, joint smoothly approaches zero but never reaches it,
    reflecting unestablished parenteral U=U.
    """
    return p_transmit(viral_load, route) * np.asarray(
        pep_efficacy_given_transmission, dtype=float
    )


def _print_table():
    routes = ['receptive_vaginal', 'receptive_anal',
              'needle_sharing', 'needlestick_hollow']
    vls = [50, 200, 1000, 10000, 30000, 100000, 1_000_000]

    header = f"{'VL (cp/mL)':>11}  " + "  ".join(f"{r:>20}" for r in routes)
    print(header)
    print("-" * len(header))
    for vl in vls:
        row = f"{vl:>11.0f}  " + "  ".join(
            f"{p_transmit(vl, r):>20.3e}" for r in routes
        )
        print(row)

    print()
    print("At VL=50:  sexual p_transmit = 0  (PARTNER2 floor enforced).")
    print(f"At VL=50:  parenteral needle_sharing p_transmit = "
          f"{p_transmit(50, 'needle_sharing'):.3e}.")
    print("This asymmetry is the publishable scientific finding.")


if __name__ == '__main__':
    print("Per-act transmission probability p_transmit(VL, route)")
    print("=" * 72)
    _print_table()
