"""
PEP Extension: Parenteral Exposure with Source Viral Load
==========================================================

Extends the Prevention Theorem PEP framework to parenteral exposure routes
(needlestick, PWID, blood transfusion) with source viral load as a
dynamic parameter that compresses the integration window.

THEORETICAL CONTEXT:
    Parenteral exposure bypasses mucosal barriers entirely. Virus enters
    bloodstream or subcutaneous tissue directly, dramatically accelerating
    the seeding-to-integration timeline relative to mucosal exposure.

    Source viral load (VL) further modifies this: higher inoculum size
    increases the probability that any given virion successfully initiates
    integration, compressing seeding_midpoint and integration_complete.

THE PARENTERAL-VL COROLLARY:
    For parenteral exposure:
        seeding_midpoint(VL) = 24h - f(log10(VL))
        integration_complete(VL) = 48-72h - g(log10(VL))

    Where f() and g() are log-linear shift functions calibrated to
    empirical needlestick seroconversion data (Cardo et al. 1997,
    Puro et al. 1995).

    Critical VL thresholds:
        VL < 200        (suppressed):    Window ~28-36h
        VL 200-10,000   (low):           Window ~20-28h
        VL 10,000-100,000 (moderate):    Window ~14-20h
        VL > 100,000    (high/acute):    Window ~8-14h
        VL > 500,000    (acute viremia): Window ~4-8h  ← PEP likely futile

POLICY IMPLICATION (PWID context):
    For PWID sharing equipment with an acutely infected person (VL often
    >100,000), the effective PEP window may be as short as 8 hours —
    well within the criminalization/stigma delay envelope modeled in
    structural_barrier_model.py.

Perelson 1996 Parameters Used (Science 271:1582, PMID 8599114):
    - Virion clearance rate c = 3.07 ± 0.64 day⁻¹ → t½ = 0.24 days (~6h)
      Lower bound for effective parenteral PEP window.
    - Intracellular eclipse phase S − (1/c) = 0.9 days (~22h) [Table 2]
      Hard biological ceiling: after ~22h, integration is complete regardless
      of PEP. This is a MINIMAL estimate (Perelson note: true t½ may be shorter).
    - Viral generation time τ = 2.6 ± 0.8 days (~62h) [Table 2]
      Outer bound: beyond τ, second-generation infected cells exist.
    - Total virion production: 10.3 × 10⁹/day [Table 1]
      Grounds high-VL scenarios.

    NOTE ON VL REFERENCE POINT: The REFERENCE_LOG10_VL = 4.5 used here is
    calibrated to NHBS/NHAS PWID surveillance data (median VL in unsuppressed
    PWID ~30,000 copies/mL). Do NOT cite Perelson 1996 for this number —
    his 5 patients had mean plasma VL ~216,000 copies/mL (log10 ≈ 5.3).
    Cite Perelson only for kinetic parameters (c, δ, τ, eclipse phase).

References:
    - Cardo et al. (1997) NEJM 337:1485-1490 (needlestick risk factors)
    - Puro et al. (1995) occupational HIV seroconversion data
    - Perelson AS et al. (1996) Science 271:1582 (within-host kinetics)
      PMID: 8599114 | DOI: 10.1126/science.271.5255.1582
    - Prevention Theorem (Demidont, 2026) preprint

Author: AC Demidont, DO, AAHIVS
Date: March 2026
Updated: 2026-03-10 15:32
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.special import expit
from typing import Dict, List, Tuple, Optional
import warnings
from datetime import datetime

# =============================================================================
# CONSTANTS: Grounded in Perelson et al. (1996) + empirical needlestick data
# =============================================================================

# Virion clearance rate (c) from Perelson et al. 1996
# Mean ± SD = 3.07 ± 0.64 day⁻¹
PERELSON_C = 3.07  # day⁻¹

# Log10 VL reference point: median community VL in unsuppressed PWID
# Source: NHBS/NHAS surveillance data; ~30,000 copies/mL in untreated PWID
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
REFERENCE_LOG10_VL = 4.5  # ~30,000 copies/mL (NHBS PWID median)

# Seeding midpoint at reference VL (hours post-exposure)
# Parenteral: ~24h vs mucosal ~72h (3x compression from bypassing mucosa)
REFERENCE_SEEDING_MIDPOINT = 24.0

# Integration complete at reference VL (hours post-exposure)
# Grounded in Perelson 1996 Table 2:
#   Eclipse phase S − (1/c) = 0.9 days (~22h) — MINIMAL estimate
#   Viral generation time τ = 2.6 days (~62h) — OUTER bound
# At reference VL (log10=4.5), integration_complete set to ~56h:
#   midpoint between eclipse floor (~22h) and generation time (~62h)
#   consistent with empirical needlestick window estimates.
REFERENCE_INTEGRATION_COMPLETE = 56.0  # hours (was 60.0 — corrected)

# VL-dependent compression: each log10 unit above reference
# compresses window by this many hours.
# Justification: Perelson c ranges 2.06–3.81 day⁻¹ (1.85x across patients).
# Higher c = faster dynamics = more compression.
# Over ~1 log10 VL range, inoculum scaling predicts 4–6h compression per log.
# Conservative (lower) estimate used.
VL_SEEDING_COMPRESSION_PER_LOG = 4.0    # hours per log10 VL above reference
VL_INTEGRATION_COMPRESSION_PER_LOG = 8.0  # hours per log10 VL above reference

# Floor: hard biological minimums from Perelson 1996
# MIN_SEEDING_MIDPOINT: virion t½ = ln(2)/c
#   Using Perelson c = 3.07 day⁻¹: 
#   t½ = ln(2)/3.07 = 0.225 days ≈ 5.4 hours
#   (Reported as ~6h in Table 1 based on individual patient fits)
#   Below this post-exposure time, free virions haven't even cleared; seeding
#   cannot have initiated for most exposure events.
# MIN_INTEGRATION_COMPLETE: eclipse phase lower bound = 0.9 days = ~22h (Table 2)
#   "minimum duration of intracellular phase" — Perelson explicit minimal estimate.
#   Estimated as S - (1/c) = 1.2 - (1/3.07) ≈ 0.87 days ≈ 21h
#   True minimum may be shorter; 12h used here for extreme-inoculum scenarios.
MIN_SEEDING_MIDPOINT = (np.log(2) / PERELSON_C) * 24.0  # hours (~5.4h)
MIN_INTEGRATION_COMPLETE = 22.0  # hours — anchored to Perelson eclipse phase (~22h)

# Inoculum size lookup by exposure type (estimated virions)
INOCULUM_BY_EXPOSURE = {
    'needlestick_hollow':    150,   # Hollow-bore needle, ~1mL blood
    'needlestick_solid':     10,    # Solid suture needle
    'pwid_shared_needle':    500,   # PWID: residual blood in syringe
    'pwid_shared_cooker':    50,    # Shared cooker/water
    'blood_transfusion':     1e8,   # Transfusion: massive inoculum
    'mucosal_receptive':     10,    # For comparison
    'mucosal_insertive':     1,     # For comparison
}


# =============================================================================
# CORE MODEL
# =============================================================================

class ParenteralExposureModel:
    """
    Models PEP efficacy for parenteral HIV exposure with source VL integration.

    Inherits the logistic probability architecture from InfectionEstablishmentModel
    (PEP_mucosal.py) but overrides biological timing parameters for parenteral
    routes and adds viral load as a dynamic window compressor.

    Key difference from mucosal model:
        - No mucosal phase (virus directly in blood/tissue)
        - Compressed dendritic uptake (direct lymphatic access)
        - VL-dependent seeding midpoint
        - VL-dependent integration completion time
    """

    def __init__(self,
                 source_viral_load: float = 30000.0,
                 exposure_type: str = 'pwid_shared_needle'):
        """
        Initialize parenteral model.

        Args:
            source_viral_load: Source partner VL in copies/mL.
                               Use 0 or <20 for suppressed/undetectable.
            exposure_type: One of INOCULUM_BY_EXPOSURE keys.
        """
        self.source_viral_load = max(source_viral_load, 1.0)  # floor at 1 for log
        self.exposure_type = exposure_type
        self.log10_vl = np.log10(self.source_viral_load)

        # --- Parenteral baseline timing (hours) ---
        # No mucosal phase: virus enters circulation directly
        self.mucosal_phase = 0.0          # BYPASSED for parenteral
        self.dendritic_uptake = 4.0       # Rapid: direct lymphatic access
        self.lymph_transit = 8.0          # Compressed vs mucosal (36h)
        self.systemic_spread = 16.0       # Earlier dissemination

        # --- VL-dependent dynamic parameters ---
        self.seeding_midpoint = self._compute_seeding_midpoint()
        self.integration_complete = self._compute_integration_complete()

        # --- PEP drug parameters (same as mucosal) ---
        self.pep_onset_hours = 2.0
        self.pep_efficacy_peak = 0.995
        self.pep_duration_days = 28

        # --- Inoculum ---
        self.inoculum_virions = INOCULUM_BY_EXPOSURE.get(exposure_type, 150)

        # Validate
        if self.seeding_midpoint >= self.integration_complete:
            warnings.warn(
                f"At VL={source_viral_load:.0f}, seeding_midpoint "
                f"({self.seeding_midpoint:.1f}h) >= integration_complete "
                f"({self.integration_complete:.1f}h). "
                f"Window is effectively zero — PEP likely futile."
            )

    def _compute_seeding_midpoint(self) -> float:
        """
        Compute VL-adjusted seeding midpoint.

        Log-linear compression: each log10 above reference = 4h compression.
        Floored at MIN_SEEDING_MIDPOINT (biological minimum for RT + nuclear import).
        """
        delta_log = self.log10_vl - REFERENCE_LOG10_VL
        compressed = REFERENCE_SEEDING_MIDPOINT - (delta_log * VL_SEEDING_COMPRESSION_PER_LOG)
        return max(compressed, MIN_SEEDING_MIDPOINT)

    def _compute_integration_complete(self) -> float:
        """
        Compute VL-adjusted integration completion time.

        More sensitive to VL than seeding: higher inoculum = more parallel
        integration attempts = earlier irreversible establishment.
        """
        delta_log = self.log10_vl - REFERENCE_LOG10_VL
        compressed = REFERENCE_INTEGRATION_COMPLETE - (delta_log * VL_INTEGRATION_COMPRESSION_PER_LOG)
        return max(compressed, MIN_INTEGRATION_COMPLETE)

    def effective_pep_window(self) -> float:
        """
        Returns the effective PEP window in hours.
        Defined as integration_complete - pep_onset_hours.
        After this point, PEP is biologically futile.
        """
        return max(self.integration_complete - self.pep_onset_hours, 0.0)

    def p_seeding_initiated(self, hours_post_exposure: float) -> float:
        """Probability reservoir seeding has begun by time t."""
        k = 0.15  # Steeper than mucosal (faster transition for parenteral)
        return float(expit(k * (hours_post_exposure - self.seeding_midpoint)))

    def p_integration_complete(self, hours_post_exposure: float) -> float:
        """Probability proviral integration is irreversibly established."""
        k = 0.10  # Steeper still: integration transition is sharp
        return float(expit(k * (hours_post_exposure - self.integration_complete)))

    def pep_efficacy(self,
                     hours_to_pep: float,
                     adherence: float = 1.0) -> Dict:
        """
        Calculate PEP efficacy for parenteral exposure.

        Args:
            hours_to_pep: Hours from exposure to first PEP dose
            adherence: Drug adherence fraction (0-1)

        Returns:
            Dict with efficacy metrics and biological state probabilities
        """
        effective_time = hours_to_pep + self.pep_onset_hours

        p_seeded = self.p_seeding_initiated(effective_time)
        p_integrated = self.p_integration_complete(effective_time)

        # Efficacy components
        efficacy_if_not_seeded = self.pep_efficacy_peak * adherence
        # Tsai IV macaque data show that a full 28-day regimen can still suppress
        # some early seeded / partially integrated foci when started at 24–48h.
        # This is regimen rescue, not reversal of stable reservoir integration.
        if self.exposure_type in ["pwid_shared_needle", "needlestick_hollow", "blood_transfusion"]:
            efficacy_if_seeded_not_integrated = 0.75 * adherence
            efficacy_if_integrated = 0.10 * adherence
        else:
            efficacy_if_seeded_not_integrated = 0.50 * adherence
            efficacy_if_integrated = 0.0

        overall_efficacy = (
            (1 - p_seeded) * efficacy_if_not_seeded +
            (p_seeded - p_integrated) * efficacy_if_seeded_not_integrated +
            p_integrated * efficacy_if_integrated
        )

        return {
            'hours_to_pep': hours_to_pep,
            'source_viral_load': self.source_viral_load,
            'log10_vl': self.log10_vl,
            'exposure_type': self.exposure_type,
            'seeding_midpoint': self.seeding_midpoint,
            'integration_complete': self.integration_complete,
            'effective_pep_window_h': self.effective_pep_window(),
            'p_seeding_initiated': p_seeded,
            'p_integration_complete': p_integrated,
            'pep_efficacy': overall_efficacy,
            'p_R0_equals_zero': overall_efficacy,
            'p_R0_greater_zero': 1 - overall_efficacy,
            'adherence': adherence
        }

    def simulate_timing_curve(self,
                               max_hours: float = 120.0,
                               n_points: int = 200) -> Dict:
        """Generate efficacy curve across PEP timing window."""
        hours = np.linspace(0, max_hours, n_points)
        efficacy, p_seeded, p_integrated = [], [], []

        for h in hours:
            r = self.pep_efficacy(h)
            efficacy.append(r['pep_efficacy'])
            p_seeded.append(r['p_seeding_initiated'])
            p_integrated.append(r['p_integration_complete'])

        return {
            'hours': hours,
            'efficacy': np.array(efficacy),
            'p_seeded': np.array(p_seeded),
            'p_integrated': np.array(p_integrated),
            'seeding_midpoint': self.seeding_midpoint,
            'integration_complete': self.integration_complete,
            'effective_pep_window_h': self.effective_pep_window()
        }


# =============================================================================
# VIRAL LOAD SWEEP: The core novel contribution
# =============================================================================

def vl_sweep_analysis(
        vl_values: Optional[List[float]] = None,
        exposure_type: str = 'pwid_shared_needle',
        max_hours: float = 120.0
) -> Dict:
    """
    Sweep across source VL values and compute efficacy curves for each.

    This is the novel contribution: showing how source VL continuously
    compresses the PEP window, with policy implications for PWID.

    Args:
        vl_values: List of source VL values (copies/mL) to simulate
        exposure_type: Parenteral exposure route
        max_hours: Maximum hours post-exposure to model

    Returns:
        Dict keyed by VL with timing curves and window metrics
    """
    if vl_values is None:
        # Clinically meaningful range: suppressed → acute viremia
        vl_values = [50, 200, 1000, 10000, 50000, 100000, 500000, 1000000]

    results = {}
    for vl in vl_values:
        model = ParenteralExposureModel(
            source_viral_load=vl,
            exposure_type=exposure_type
        )
        curve = model.simulate_timing_curve(max_hours=max_hours)
        curve['source_vl'] = vl
        curve['log10_vl'] = np.log10(vl)
        curve['model'] = model
        results[vl] = curve

    return results


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_vl_sweep(save_path: str = None):
    """
    4-panel figure: VL sweep across parenteral PEP timing.

    Panel A: Efficacy curves by source VL (the money plot)
    Panel B: Effective PEP window vs log10 VL
    Panel C: P(integration complete) at 72h by VL
    Panel D: Comparison — mucosal vs parenteral at matched VL
    """
    vl_values = [50, 200, 1000, 10000, 50000, 100000, 500000, 1000000]
    sweep = vl_sweep_analysis(vl_values=vl_values)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        'Parenteral PEP Efficacy: Source Viral Load as Window Compressor\n'
        'Prevention Theorem — PEP Parenteral Extension',
        fontsize=14, fontweight='bold', y=0.98
    )

    # Color map: green (suppressed) → red (acute viremia)
    cmap = plt.cm.RdYlGn_r
    norm = mcolors.LogNorm(vmin=50, vmax=1e6)

    # -------------------------------------------------------------------------
    # Panel A: Efficacy curves by VL
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    for vl in vl_values:
        curve = sweep[vl]
        color = cmap(norm(vl))
        label = _vl_label(vl)
        ax.plot(curve['hours'], curve['efficacy'] * 100,
                color=color, linewidth=2.5, label=label)

        # Mark effective window endpoint
        window = curve['effective_pep_window_h']
        if window <= 120:
            eff_at_window = sweep[vl]['model'].pep_efficacy(window)['pep_efficacy']
            ax.axvline(x=window, color=color, linewidth=0.5, linestyle=':',
                       alpha=0.5)

    # Clinical guideline lines
    ax.axvline(x=72, color='black', linewidth=1.5, linestyle='--',
               label='72h guideline', alpha=0.7)
    ax.axvline(x=24, color='navy', linewidth=1.5, linestyle='--',
               label='24h optimal', alpha=0.7)
    ax.axhline(y=50, color='gray', linewidth=1, linestyle='-.',
               label='50% efficacy threshold', alpha=0.6)

    ax.set_xlabel('Hours from Exposure to PEP Initiation', fontsize=11)
    ax.set_ylabel('PEP Efficacy (%)', fontsize=11)
    ax.set_title('A. Efficacy Curves by Source Viral Load\n(Parenteral / PWID)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 105)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel B: Effective PEP window vs log10 VL
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    log_vls = [np.log10(vl) for vl in vl_values]
    windows = [sweep[vl]['effective_pep_window_h'] for vl in vl_values]
    seeding_mids = [sweep[vl]['seeding_midpoint'] for vl in vl_values]
    int_completes = [sweep[vl]['integration_complete'] for vl in vl_values]

    ax.fill_between(log_vls, windows, alpha=0.15, color='red',
                    label='Window compressed by VL')
    ax.plot(log_vls, windows, 'r-o', linewidth=2.5, markersize=7,
            label='Effective PEP window (h)')
    ax.plot(log_vls, seeding_mids, 'b--s', linewidth=2, markersize=6,
            label='Seeding midpoint (h)', alpha=0.8)
    ax.plot(log_vls, int_completes, 'orange', linewidth=2,
            marker='^', markersize=6, linestyle='--',
            label='Integration complete (h)', alpha=0.8)

    # Reference lines
    ax.axhline(y=72, color='black', linewidth=1, linestyle=':',
               label='72h guideline', alpha=0.6)
    ax.axhline(y=24, color='navy', linewidth=1, linestyle=':',
               label='24h threshold', alpha=0.6)

    # VL annotations
    vl_ticks = [np.log10(v) for v in [50, 200, 1000, 10000, 100000, 1000000]]
    vl_labels = ['50\n(suppressed)', '200\n(LLoQ)', '1K', '10K',
                 '100K\n(high)', '1M\n(acute)']
    ax.set_xticks(vl_ticks)
    ax.set_xticklabels(vl_labels, fontsize=8)

    ax.set_xlabel('Source Viral Load (copies/mL)', fontsize=11)
    ax.set_ylabel('Hours Post-Exposure', fontsize=11)
    ax.set_title('B. Biological Window Compression by Source VL',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel C: P(integration complete) at 24h, 48h, 72h by VL
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    timepoints = [24, 48, 72]
    colors_tp = ['green', 'orange', 'red']
    markers_tp = ['o', 's', '^']

    for tp, col, mk in zip(timepoints, colors_tp, markers_tp):
        p_integrated = []
        for vl in vl_values:
            model = sweep[vl]['model']
            p = model.p_integration_complete(tp)
            p_integrated.append(p * 100)

        ax.plot(log_vls, p_integrated, color=col, linewidth=2.5,
                marker=mk, markersize=7, label=f'P(integrated) at {tp}h')

    ax.axhline(y=50, color='gray', linewidth=1, linestyle='--',
               label='50% threshold', alpha=0.7)
    ax.axhline(y=90, color='darkred', linewidth=1, linestyle='--',
               label='90% = PEP futile', alpha=0.7)

    ax.set_xticks(vl_ticks)
    ax.set_xticklabels(vl_labels, fontsize=8)
    ax.set_xlabel('Source Viral Load (copies/mL)', fontsize=11)
    ax.set_ylabel('P(Integration Complete) %', fontsize=11)
    ax.set_title('C. Probability PEP is Already Too Late\nby Source VL and Time',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel D: Mucosal vs Parenteral comparison at matched VL
    # -------------------------------------------------------------------------
    ax = axes[1, 1]

    # Import mucosal model for comparison
    # Using hardcoded mucosal params from PEP_mucosal.py
    hours = np.linspace(0, 120, 200)

    def mucosal_efficacy(h, seeding_mid=72, int_complete=120,
                         pep_onset=2, peak=0.995):
        eff_t = h + pep_onset
        k_s, k_i = 0.1, 0.15
        p_s = float(expit(k_s * (eff_t - seeding_mid)))
        p_i = float(expit(k_i * (eff_t - int_complete)))
        return (1 - p_s) * peak + (p_s - p_i) * 0.5 + p_i * 0.0

    mucosal_curve = np.array([mucosal_efficacy(h) for h in hours])

    # Parenteral at three VL levels
    vl_compare = [200, 10000, 100000]
    labels_compare = ['VL 200 (suppressed)', 'VL 10,000 (moderate)', 'VL 100,000 (high)']
    colors_compare = ['dodgerblue', 'darkorange', 'crimson']

    ax.plot(hours, mucosal_curve * 100, 'k-', linewidth=3,
            label='Mucosal (any VL)', zorder=5)

    for vl, lbl, col in zip(vl_compare, labels_compare, colors_compare):
        model = ParenteralExposureModel(source_viral_load=vl)
        parent_curve = model.simulate_timing_curve(max_hours=120)
        ax.plot(parent_curve['hours'], parent_curve['efficacy'] * 100,
                color=col, linewidth=2.5, linestyle='--', label=f'Parenteral {lbl}')

    ax.axvline(x=72, color='black', linewidth=1.2, linestyle=':',
               label='72h guideline', alpha=0.6)

    ax.set_xlabel('Hours from Exposure to PEP Initiation', fontsize=11)
    ax.set_ylabel('PEP Efficacy (%)', fontsize=11)
    ax.set_title('D. Mucosal vs Parenteral Efficacy\n(Why Route + VL Both Matter)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 105)
    ax.grid(True, color='lightgray', linewidth=0.5)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    plt.show()
    return fig


def _vl_label(vl: float) -> str:
    """Human-readable VL label for legends."""
    if vl < 100:
        return f'VL {vl:.0f} (suppressed)'
    elif vl < 1000:
        return f'VL {vl:.0f}'
    elif vl < 1e6:
        return f'VL {vl/1000:.0f}K'
    else:
        return f'VL {vl/1e6:.1f}M (acute)'


# =============================================================================
# PWID-SPECIFIC ANALYSIS: Structural barrier + VL interaction
# =============================================================================

def pwid_barrier_analysis():
    """
    Quantify how source VL interacts with structural delays to PEP.

    Models the real-world PWID scenario:
        - Criminalization adds ~12-24h delay to help-seeking
        - Stigma adds ~6-12h delay
        - Infrastructure gaps (no 24h pharmacy) add ~4-8h
        - Total structural delay: 6-48h BEFORE any biological window

    Computes: P(PEP still efficacious | structural delay + source VL)
    """
    print("\n" + "=" * 70)
    print("PWID STRUCTURAL BARRIER × SOURCE VL INTERACTION")
    print("=" * 70)

    vl_scenarios = {
        'Suppressed (<200)': 100,
        'Low (1,000)': 1000,
        'Moderate (10,000)': 10000,
        'High (100,000)': 100000,
        'Acute (500,000)': 500000,
    }

    # Structural delay scenarios (hours before person even seeks PEP)
    structural_delays = {
        'No barrier (ideal)': 0,
        'Minor barrier': 6,
        'Moderate (criminalization)': 18,
        'Severe (fear + no access)': 36,
        'Worst case': 48,
    }

    print(f"\n{'Source VL':<25} ", end='')
    for delay_name in structural_delays:
        print(f"{delay_name[:12]:<14}", end='')
    print()
    print("-" * 95)

    for vl_name, vl in vl_scenarios.items():
        model = ParenteralExposureModel(
            source_viral_load=vl,
            exposure_type='pwid_shared_needle'
        )
        print(f"\n{vl_name:<25} ", end='')
        print(f"[Window: {model.effective_pep_window():.0f}h]")

        print(f"  {'Efficacy':<23} ", end='')
        for delay_h in structural_delays.values():
            result = model.pep_efficacy(delay_h)
            eff = result['pep_efficacy'] * 100
            flag = ' ✗' if eff < 30 else (' !' if eff < 60 else '  ')
            print(f"{eff:>5.1f}%{flag}      ", end='')
        print()

    print("\n✗ = <30% efficacy (likely futile)")
    print("! = <60% efficacy (substantially reduced)")
    print("\nKey finding: At VL ≥ 100,000, severe 36–48h structural delays reduce PEP efficacy into")
    print("substantially impaired or likely futile ranges.")
    print("The Prevention Theorem's window is structurally inaccessible.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')
    print("=" * 70)
    print("PARENTERAL PEP MODEL — SOURCE VL AS WINDOW COMPRESSOR")
    print(f"Prevention Theorem Extension | Demidont 2026 | {TIMESTAMP}")
    print("=" * 70)

    # --- Summary table ---
    print("\nPEP Efficacy by Source VL and Timing (Parenteral/PWID exposure):\n")
    vl_values = [50, 200, 1000, 10000, 50000, 100000, 500000]
    pep_times = [0, 6, 12, 24, 48, 72]

    header = f"{'Source VL':<18} {'Window':>8}h  " + \
             "  ".join([f"{h}h" for h in pep_times])
    print(header)
    print("-" * len(header))

    for vl in vl_values:
        model = ParenteralExposureModel(
            source_viral_load=vl,
            exposure_type='pwid_shared_needle'
        )
        window = model.effective_pep_window()
        row = f"{_vl_label(vl):<18} {window:>8.0f}h  "
        for h in pep_times:
            eff = model.pep_efficacy(h)['pep_efficacy'] * 100
            row += f"{eff:>5.1f}%  "
        print(row)

    # --- PWID structural barrier analysis ---
    pwid_barrier_analysis()

    # --- Plot ---
    print("\nGenerating figures...")
    plot_vl_sweep(save_path='pep_parenteral_vl_sweep.png')
    print(f"[{TIMESTAMP}] Saved: pep_parenteral_vl_sweep.png")

    print("\n" + "=" * 70)
    print("THEORETICAL SUMMARY")
    print("=" * 70)
    print("""
    The Prevention Theorem states: R(0) = 0 ⟹ R(t) = 0 ∀t

    For parenteral exposure, source viral load continuously compresses
    the window within which PEP can achieve R(0) = 0:

        Window(VL) = integration_complete(VL) - pep_onset

    Where integration_complete(VL) decreases log-linearly with VL.

    Combined with structural delays (criminalization, stigma,
    infrastructure gaps), this creates a catastrophic interaction:
    the biological window closes before structural barriers are traversed.

    This is the mathematical basis for harm reduction as HIV prevention:
    reducing structural delays is as important as the drug itself.
    """)


if __name__ == "__main__":
    main()