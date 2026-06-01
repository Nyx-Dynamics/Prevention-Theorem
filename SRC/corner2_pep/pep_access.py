"""
P_access — Barrier-Gated PEP Access Probability
=================================================

Computes the probability that an individual with a given structural
barrier profile successfully reaches PEP care within the finite
kinetic window (i.e., before eclipse_boundary_hours).

Design:
  - Mirrors the BTG adjusted_success_rate calculation but applied to
    PEP access, not PrEP bridge-period retention.
  - Uses only access-relevant barriers (TRANSPORTATION, INSURANCE_DELAYS,
    etc.) — NOT long-term adherence barriers.
  - Eclipse boundary is a hard floor; P_access × efficacy(t) → 0 beyond it.
  - Double-counting guard: P_access operates on exposures NOT already
    averted by PrEP (averted_by_prep flag or explicit complement).

Import kinetic constants ONLY from substrate_fpw — never re-define them.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from corner2_pep.substrate_fpw import (
    FPW_SUBSTRATE,
    PEP_ACCESS_WEIGHTS,
    FPWKineticSubstrate,
    PEPAccessBarrierWeights,
)


# ---------------------------------------------------------------------------
# BARRIER PROFILE
# ---------------------------------------------------------------------------

@dataclass
class BarrierProfile:
    """
    Active barriers for a patient / population stratum.

    barrier_keys: list of BTG UPPERCASE_SNAKE_CASE barrier keys present
    averted_by_prep: fraction of exposures already averted by PrEP success.
                     P_access only applies to the complement (1 - averted_by_prep).
    """
    barrier_keys: List[str] = field(default_factory=list)
    averted_by_prep: float = 0.0      # double-counting guard
    population_key: str = "GENERAL"


# ---------------------------------------------------------------------------
# P_ACCESS CALCULATION
# ---------------------------------------------------------------------------

def p_access(
    profile: BarrierProfile,
    substrate: FPWKineticSubstrate = FPW_SUBSTRATE,
    weights: PEPAccessBarrierWeights = PEP_ACCESS_WEIGHTS,
    diminishing_returns_factor: float = 0.70,
    max_attrition_ceiling: float = 0.95,
) -> Dict:
    """
    Compute P_access: probability of reaching PEP care within the window
    given the barrier profile.

    Parameters
    ----------
    profile : BarrierProfile
    substrate : FPWKineticSubstrate — kinetics substrate (immutable)
    weights : PEPAccessBarrierWeights — barrier weight lookup
    diminishing_returns_factor : float — per-barrier stacking diminishing returns
    max_attrition_ceiling : float — maximum achievable access-attrition

    Returns
    -------
    dict with keys:
      p_access_raw        : P(reach PEP in window) ignoring averted exposures
      p_access_net        : corrected for averted_by_prep (double-counting guard)
      access_attrition    : total access barrier impact
      active_barriers     : subset of profile.barrier_keys that are access-relevant
      inactive_barriers   : profile barriers not in the access-barrier set
    """
    access_barrier_set = set(weights.ACCESS_BARRIER_KEYS)
    active = [k for k in profile.barrier_keys if k in access_barrier_set]
    inactive = [k for k in profile.barrier_keys if k not in access_barrier_set]

    # Baseline access probability = 1.0 (no barriers → always reaches PEP in window)
    # Attrition accumulates with each active barrier, with diminishing returns
    access_attrition = 0.0
    for i, key in enumerate(active):
        w = getattr(weights, key, 0.0)
        # Diminishing returns: each successive barrier has less marginal impact
        factor = diminishing_returns_factor ** i
        access_attrition += w * factor

    access_attrition = min(access_attrition, max_attrition_ceiling)
    p_raw = max(0.0, 1.0 - access_attrition)

    # Double-counting guard: PEP only recovers exposures NOT already averted by PrEP
    # averted_by_prep = adjusted_success_rate from the BTG model for this population
    exposures_not_averted = max(0.0, 1.0 - profile.averted_by_prep)
    p_net = p_raw * exposures_not_averted

    return {
        "p_access_raw": p_raw,
        "p_access_net": p_net,
        "access_attrition": access_attrition,
        "active_barriers": active,
        "inactive_barriers": inactive,
        "averted_by_prep": profile.averted_by_prep,
        "exposures_not_averted": exposures_not_averted,
        "population_key": profile.population_key,
        "barrier_source": "BTG v2.2.0 PRE-WISE",
    }
