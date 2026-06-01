"""
FPW Kinetic Substrate — Finite Prevention Window Corner 2
==========================================================

SINGLE SOURCE OF TRUTH for all kinetic constants used in the PEP
finite-window layer. No other module may define these constants directly;
the pre-commit hook enforces this. Import from here.

References (do not cite Perelson for VL reference point — see notes):
  Perelson AS et al. Science 1996;271:1582. PMID 8599114.
    → kinetic parameters: c, delta, tau_eclipse, viral generation time
  Cardo et al. NEJM 1997;337:1485.
    → seeding/integration calibration for parenteral exposure
  NHBS/NHAS PWID surveillance
    → REFERENCE_LOG10_VL = 4.5 (~30,000 copies/mL median unsuppressed PWID)
      NOT from Perelson (his 5 patients: mean VL ~216,000 = log10 5.3)

Provenance: Prevention-Theorem repo. Constants previously hardcoded in
PEP_parenteral_perelson.py module scope; centralized here for Corner-2.
"""

from __future__ import annotations
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# SUBSTRATE: FPW KINETICS
# These are the ONLY authorised definitions of these constants.
# Pre-commit hook blocks re-definition outside this file.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FPWKineticSubstrate:
    """
    Immutable kinetic substrate for the Finite Prevention Window model.

    All values are Perelson-1996-anchored (parenteral/IV route baseline).
    VL-dependent compression is handled by the ParenteralExposureModel
    in perelson/parenteral/PEP_parenteral_perelson.py using these as base.
    """
    # --- Eclipse phase (Perelson 1996, Table 2) ---
    tau_eclipse_days: float = 0.9          # ~21.6 h; biological ceiling for PEP
    tau_eclipse_hours: float = 21.6        # convenience alias

    # --- Virion clearance (Perelson 1996, Table 2) ---
    c_per_day: float = 23.0                # clearance rate (day⁻¹)
    virion_half_life_hours: float = 0.72   # ln(2)/c ≈ 0.72 h

    # --- Integration window (parenteral, reference VL) ---
    # At REFERENCE_LOG10_VL = 4.5 (~30,000 copies/mL, NHBS PWID median)
    reference_log10_vl: float = 4.5        # NHBS/NHAS PWID surveillance; NOT Perelson
    reference_seeding_midpoint_hours: float = 24.0
    reference_integration_complete_hours: float = 56.0  # corrected from 60.0

    # --- PEP efficacy ceiling ---
    pep_efficacy_peak: float = 0.995       # maximum biological efficacy (perfect adherence, t→0)

    # --- VL compression coefficients ---
    vl_seeding_compression_per_log: float = 4.0    # hours compressed per log10 VL unit
    vl_integration_compression_per_log: float = 8.0

    # --- Operational floor: eclipse boundary enforcement ---
    # PEP efficacy is set to 0 at or beyond tau_eclipse_hours.
    # This is the hard biological ceiling; no access factor can recover it.
    eclipse_boundary_hours: float = 21.6   # == tau_eclipse_hours; explicit for clarity

    # --- Viral generation time (Perelson 1996, Table 2) ---
    viral_generation_time_days: float = 2.6   # ~62 h; outer bound for second-generation cells


# Singleton: import this rather than constructing per-call
FPW_SUBSTRATE = FPWKineticSubstrate()


@dataclass(frozen=True)
class PEPAccessBarrierWeights:
    """
    Access-relevant BTG barrier weights for PEP gate calculation.

    Source: lai-prep-bridge-tool-pub SCR/code/algorithm/lai_prep_config.json
    config v2.2.0. Only barriers that specifically impede reaching PEP care
    within the finite kinetic window are included — long-term adherence
    barriers (CHILDCARE, COMPETING_PRIORITIES) are excluded.

    Label: PRE-WISE — these weights are literature-grounded (BTG v2.2).
    WISE-fitted weights will supersede this once WISE data are available.
    """
    TRANSPORTATION: float = 0.10
    INSURANCE_DELAYS: float = 0.12
    SCHEDULING_CONFLICTS: float = 0.05
    MEDICAL_MISTRUST: float = 0.10
    HEALTHCARE_DISCRIMINATION: float = 0.12
    LEGAL_CONCERNS: float = 0.15
    LACK_IDENTIFICATION: float = 0.10
    PRIVACY_CONCERNS: float = 0.08

    # Access-barrier keys in evaluation order (most time-sensitive first)
    ACCESS_BARRIER_KEYS: tuple = field(default_factory=lambda: (
        "TRANSPORTATION",
        "INSURANCE_DELAYS",
        "SCHEDULING_CONFLICTS",
        "MEDICAL_MISTRUST",
        "HEALTHCARE_DISCRIMINATION",
        "LEGAL_CONCERNS",
        "LACK_IDENTIFICATION",
        "PRIVACY_CONCERNS",
    ))


PEP_ACCESS_WEIGHTS = PEPAccessBarrierWeights()
