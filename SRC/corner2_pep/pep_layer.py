"""
PEP Finite-Window Layer — In-Window Recovery, Net Protection, Deficit
======================================================================

Implements the second corner of the Prevention Theorem:

    in_window_recovery(t, profile) = P_access(profile) × efficacy(t)

where efficacy(t) declines to 0 at the eclipse boundary.

Net additive protection:
    net_protection = PrEP_adjusted_success + in_window_recovery
        (capped at 1.0; double-counting guard already in P_access)

Deficit:
    potential = P_access_raw × efficacy_at_best_access (t → 0)
    realized  = in_window_recovery
    deficit   = potential - realized

All kinetic constants imported from substrate_fpw. Never re-defined here.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from typing import Dict, List, Optional

# Path resolution — allow import whether run from repo root or SRC/
_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from corner2_pep.substrate_fpw import FPW_SUBSTRATE, FPWKineticSubstrate
from corner2_pep.pep_access import BarrierProfile, p_access, PEP_ACCESS_WEIGHTS

# PEP kinetics model from the existing parenteral module
from perelson.parenteral.PEP_parenteral_perelson import ParenteralExposureModel


# ---------------------------------------------------------------------------
# IN-WINDOW EFFICACY (wraps existing ParenteralExposureModel)
# ---------------------------------------------------------------------------

def pep_efficacy_at_delay(
    hours_to_pep: float,
    source_log10_vl: float = None,
    substrate: FPWKineticSubstrate = FPW_SUBSTRATE,
) -> float:
    """
    PEP biological efficacy at a given delay post-exposure.

    Eclipse boundary is enforced: efficacy = 0 at or beyond
    substrate.eclipse_boundary_hours (hard biological ceiling).

    Parameters
    ----------
    hours_to_pep : float — hours from exposure to PEP initiation
    source_log10_vl : float — source log10 viral load; defaults to
                      substrate.reference_log10_vl if None
    substrate : FPWKineticSubstrate

    Returns
    -------
    float in [0, 1]
    """
    # Hard eclipse boundary
    if hours_to_pep >= substrate.eclipse_boundary_hours:
        return 0.0

    vl = source_log10_vl if source_log10_vl is not None else substrate.reference_log10_vl
    source_vl_copies = 10.0 ** vl
    model = ParenteralExposureModel(source_viral_load=source_vl_copies)
    result = model.pep_efficacy(hours_to_pep)
    return float(result.get("pep_efficacy", 0.0))


# ---------------------------------------------------------------------------
# IN-WINDOW RECOVERY
# ---------------------------------------------------------------------------

def in_window_recovery(
    profile: BarrierProfile,
    hours_to_pep: float,
    source_log10_vl: float = None,
    substrate: FPWKineticSubstrate = FPW_SUBSTRATE,
) -> Dict:
    """
    Compute in-window recovery for a given barrier profile and PEP delay.

    in_window_recovery = P_access(profile) × efficacy(hours_to_pep)

    The double-counting guard (averted_by_prep) is already applied inside
    p_access — no additional correction needed here.

    Returns
    -------
    dict with full breakdown including P_access, efficacy, recovery, deficit
    """
    access = p_access(profile, substrate=substrate)
    eff = pep_efficacy_at_delay(hours_to_pep, source_log10_vl, substrate)

    recovery_raw = access["p_access_raw"] * eff
    recovery_net = access["p_access_net"] * eff

    # Potential = what recovery would be with zero barriers at this delay
    # (P_access_raw = 1.0, same efficacy)
    potential = eff  # no barrier attrition
    deficit = max(0.0, potential - recovery_raw)

    return {
        # --- Inputs ---
        "hours_to_pep": hours_to_pep,
        "source_log10_vl": source_log10_vl or substrate.reference_log10_vl,
        "eclipse_boundary_hours": substrate.eclipse_boundary_hours,
        "eclipse_exceeded": hours_to_pep >= substrate.eclipse_boundary_hours,
        # --- P_access breakdown ---
        **{f"access_{k}": v for k, v in access.items()},
        # --- Efficacy ---
        "pep_efficacy_at_delay": eff,
        # --- Recovery ---
        "in_window_recovery_raw": recovery_raw,   # ignoring PrEP-averted
        "in_window_recovery_net": recovery_net,   # net of PrEP-averted
        # --- Deficit ---
        "potential_recovery": potential,
        "deficit": deficit,
        "deficit_pct_of_potential": (deficit / potential * 100) if potential > 0 else 0.0,
    }


# ---------------------------------------------------------------------------
# ADDITIVE NET PROTECTION
# ---------------------------------------------------------------------------

def net_protection(
    prep_adjusted_success: float,
    profile: BarrierProfile,
    hours_to_pep: float,
    source_log10_vl: float = None,
    substrate: FPWKineticSubstrate = FPW_SUBSTRATE,
) -> Dict:
    """
    Additive net protection from both PrEP persistence and PEP in-window recovery.

    Net = PrEP_adjusted_success + in_window_recovery_net
    Capped at 1.0 (cannot exceed full protection).

    The double-counting guard is embedded: profile.averted_by_prep should be
    set to prep_adjusted_success so P_access only applies to the complement.

    Parameters
    ----------
    prep_adjusted_success : float — adjusted_success_rate from BTG model
    profile : BarrierProfile — must have averted_by_prep == prep_adjusted_success
    hours_to_pep : float
    source_log10_vl : float
    substrate : FPWKineticSubstrate
    """
    assert abs(profile.averted_by_prep - prep_adjusted_success) < 1e-9, (
        "double-counting guard: profile.averted_by_prep must equal "
        "prep_adjusted_success to avoid counting PrEP-averted exposures twice"
    )

    recovery = in_window_recovery(profile, hours_to_pep, source_log10_vl, substrate)
    combined = min(1.0, prep_adjusted_success + recovery["in_window_recovery_net"])

    return {
        "prep_adjusted_success": prep_adjusted_success,
        "pep_in_window_recovery_net": recovery["in_window_recovery_net"],
        "net_protection": combined,
        "net_protection_gain_over_prep_alone": combined - prep_adjusted_success,
        "pep_deficit": recovery["deficit"],
        "pep_deficit_pct_of_potential": recovery["deficit_pct_of_potential"],
        "population_key": profile.population_key,
        "hours_to_pep": hours_to_pep,
        "eclipse_exceeded": recovery["eclipse_exceeded"],
        "barrier_source": recovery.get("access_barrier_source", "BTG v2.2.0 PRE-WISE"),
        # P_access breakdown (pass-through for run scripts)
        "p_access_raw": recovery["access_p_access_raw"],
        "p_access_net": recovery["access_p_access_net"],
        "access_attrition": recovery["access_access_attrition"],
        "active_barriers": recovery["access_active_barriers"],
        "pep_efficacy_at_delay": recovery["pep_efficacy_at_delay"],
    }
