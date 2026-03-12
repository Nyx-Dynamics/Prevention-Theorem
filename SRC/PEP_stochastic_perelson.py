"""
PEP Stochastic Extension: VL Uncertainty + Population-Level Efficacy
=====================================================================

Adds a stochastic layer to PEP_parenteral_perelson.py addressing the core
clinical reality: at the moment of exposure, source VL is UNKNOWN.

The clinician prescribing PEP, the PWID deciding whether to seek care,
and the public health model predicting outcomes — all operate under
uncertainty about the source's viral load.

This module computes:
    1. Expected PEP efficacy under VL uncertainty (Monte Carlo)
    2. Credible intervals on efficacy by delay time
    3. The "VL knowledge premium" — how much efficacy is gained
       if source VL is known (e.g., from a rapid test at PEP initiation)

THEORETICAL SIGNIFICANCE:
    The deterministic model (PEP_parenteral_perelson.py) gives:
        Efficacy(t, VL) — exact, given known VL

    This module gives:
        E[Efficacy(t)] = integral[ Efficacy(t, VL) * P(VL) ] dVL

    Where P(VL) is the empirical community VL distribution.

    The gap between E[Efficacy | VL known] vs E[Efficacy | VL unknown]
    is the information value of rapid source VL testing at PEP initiation.
    This has direct policy implications for ED-based PEP programs.

VL DISTRIBUTION PARAMETERS (empirical basis):
    Community VL distribution among untreated PWID approximates log-normal.
    Source: NHBS/NHAS surveillance data; Metzger et al.; Degenhardt et al.
        Mean log10 VL ~4.3–4.5 (unsuppressed PWID)
        SD log10 VL ~1.1 (wide: ranges from suppressed to >6 log)

    KINETIC PARAMETERS (Perelson et al., Science 1996, PMID 8599114):
        Virion t½ = 0.24 days (~6h)  [Table 1: c = 3.07 day⁻¹]
        Eclipse phase = 0.9 days (~22h)  [Table 2: S − 1/c, minimal estimate]
        Generation time τ = 2.6 days (~62h)  [Table 2]
    These kinetic constants ground the integration timeline in
    PEP_parenteral_perelson.py — NOT the VL distribution center.
    DO NOT conflate: cite Perelson for kinetics, NHBS for VL distribution.

Author: AC Demidont, DO, AAHIVS
Date: March 2026
Updated: 2026-03-10 15:32
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime
from scipy import stats
from scipy.special import expit
from typing import Dict, List, Optional
import warnings

# Import the deterministic model
from PEP_parenteral_perelson import ParenteralExposureModel, REFERENCE_LOG10_VL

# =============================================================================
# VL DISTRIBUTION PARAMETERS
# =============================================================================

VL_DIST_PARAMS = {
    # mean_log10 = 4.5 → ~30,000 copies/mL
    # Source: NHBS PWID surveillance median for unsuppressed individuals.
    # Perelson 1996 patients had mean ~216,000 (log10≈5.3) — cite for kinetics only.
    'pwid_untreated': {
        'mean_log10': 4.5,
        'sd_log10': 1.1,
        'label': 'PWID (untreated/unsuppressed)',
        'color': 'crimson'
    },
    'pwid_mixed_treatment': {
        'mean_log10': 3.2,
        'sd_log10': 1.5,
        'label': 'PWID (mixed treatment status)',
        'color': 'darkorange'
    },
    'general_community': {
        'mean_log10': 3.8,
        'sd_log10': 1.2,
        'label': 'General HIV+ community',
        'color': 'steelblue'
    },
    # mean_log10 = 5.0 → ~100,000 copies/mL
    # Acute HIV: typical peak viremia 10⁵–10⁷ copies/mL; Perelson patient 105
    # had VL 643,000 (log10≈5.8), consistent with this range.
    'acute_infection_enriched': {
        'mean_log10': 5.0,
        'sd_log10': 0.9,
        'label': 'Acute-infection-enriched network',
        'color': 'darkred'
    }
}


# =============================================================================
# MONTE CARLO ENGINE
# =============================================================================

class StochasticPEPModel:
    """
    Monte Carlo simulation of PEP efficacy under source VL uncertainty.

    For each PEP delay time, samples n_simulations VL values from the
    community distribution and computes the resulting efficacy distribution.
    """

    def __init__(self,
                 vl_distribution: str = 'pwid_untreated',
                 exposure_type: str = 'pwid_shared_needle',
                 n_simulations: int = 10000,
                 random_seed: int = 42):
        self.vl_dist_params = VL_DIST_PARAMS[vl_distribution]
        self.vl_distribution = vl_distribution
        self.exposure_type = exposure_type
        self.n_simulations = n_simulations
        self.rng = np.random.default_rng(random_seed)

        # Pre-sample VL values once — reuse across all time points
        self.sampled_log10_vl = self.rng.normal(
            loc=self.vl_dist_params['mean_log10'],
            scale=self.vl_dist_params['sd_log10'],
            size=n_simulations
        )
        # Clip: floor at log10(20)=1.3 (below LLoQ = effectively suppressed)
        self.sampled_log10_vl = np.clip(self.sampled_log10_vl, 1.3, 7.0)
        self.sampled_vl = 10 ** self.sampled_log10_vl

    def simulate_at_delay(self, hours_to_pep: float) -> Dict:
        """Run Monte Carlo at a single delay time."""
        efficacies = np.zeros(self.n_simulations)

        for i, vl in enumerate(self.sampled_vl):
            model = ParenteralExposureModel(
                source_viral_load=vl,
                exposure_type=self.exposure_type
            )
            efficacies[i] = model.pep_efficacy(hours_to_pep)['pep_efficacy']

        return {
            'hours_to_pep': hours_to_pep,
            'efficacies': efficacies,
            'mean': np.mean(efficacies),
            'median': np.median(efficacies),
            'ci_95_lower': np.percentile(efficacies, 2.5),
            'ci_95_upper': np.percentile(efficacies, 97.5),
            'ci_50_lower': np.percentile(efficacies, 25),
            'ci_50_upper': np.percentile(efficacies, 75),
            'p_below_50pct': np.mean(efficacies < 0.5),
            'p_below_30pct': np.mean(efficacies < 0.3),
            'p_futile': np.mean(efficacies < 0.1),
        }

    def simulate_timing_curve(self,
                               max_hours: float = 120.0,
                               n_timepoints: int = 50) -> Dict:
        """Run Monte Carlo across full timing curve."""
        hours = np.linspace(0, max_hours, n_timepoints)
        results = {k: [] for k in ['hours', 'mean', 'median',
                                    'ci_95_lower', 'ci_95_upper',
                                    'ci_50_lower', 'ci_50_upper',
                                    'p_below_50pct', 'p_below_30pct',
                                    'p_futile']}
        results['hours'] = hours

        label = self.vl_dist_params['label']
        print(f"  Running {self.n_simulations:,} sims x "
              f"{n_timepoints} timepoints ({label})...")

        for h in hours:
            r = self.simulate_at_delay(h)
            for key in ['mean', 'median', 'ci_95_lower', 'ci_95_upper',
                        'ci_50_lower', 'ci_50_upper', 'p_below_50pct',
                        'p_below_30pct', 'p_futile']:
                results[key].append(r[key])

        for key in results:
            if key != 'hours':
                results[key] = np.array(results[key])

        return results


# =============================================================================
# VL KNOWLEDGE PREMIUM
# =============================================================================

def compute_vl_knowledge_premium(
        hours_to_pep_values: List[float] = None,
        vl_distribution: str = 'pwid_untreated',
        n_simulations: int = 10000
) -> List[Dict]:
    """
    Compute the 'VL knowledge premium' at each delay time.

    Premium(t) = E[Efficacy | VL known, urgency acted upon]
               - E[Efficacy | VL unknown]

    'Known and acted upon': rapid VL test showing high VL triggers
    2h urgency reduction in effective delay for high-VL cases (VL > 10K).

    Quantifies the value of point-of-care VL testing at PEP initiation.
    """
    if hours_to_pep_values is None:
        hours_to_pep_values = [0, 6, 12, 18, 24, 36, 48, 72]

    dist_params = VL_DIST_PARAMS[vl_distribution]
    rng = np.random.default_rng(42)

    sampled_log10_vl = rng.normal(
        loc=dist_params['mean_log10'],
        scale=dist_params['sd_log10'],
        size=n_simulations
    )
    sampled_log10_vl = np.clip(sampled_log10_vl, 1.3, 7.0)
    sampled_vl = 10 ** sampled_log10_vl

    results = []
    for h in hours_to_pep_values:
        efficacies_unknown, efficacies_known = [], []

        for vl in sampled_vl:
            # Unknown VL scenario: standard delay
            eff_unknown = ParenteralExposureModel(
                source_viral_load=vl
            ).pep_efficacy(h)['pep_efficacy']
            efficacies_unknown.append(eff_unknown)

            # Known VL scenario: 2h urgency reduction if VL > 10K
            urgency_reduction = 2.0 if vl > 10000 else 0.0
            effective_h = max(h - urgency_reduction, 0)
            eff_known = ParenteralExposureModel(
                source_viral_load=vl
            ).pep_efficacy(effective_h)['pep_efficacy']
            efficacies_known.append(eff_known)

        mean_unknown = np.mean(efficacies_unknown)
        mean_known = np.mean(efficacies_known)
        premium = mean_known - mean_unknown

        results.append({
            'hours': h,
            'mean_efficacy_unknown_vl': mean_unknown,
            'mean_efficacy_known_vl': mean_known,
            'premium_absolute': premium,
            'premium_relative': premium / mean_unknown if mean_unknown > 0 else 0,
        })

    return results


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_stochastic_analysis(save_path: str = None):
    """
    4-panel stochastic figure.

    Panel A: Mean + 95% CI efficacy curves by VL distribution
    Panel B: VL distribution densities
    Panel C: P(efficacy < 50%) by delay — clinical failure probability
    Panel D: VL knowledge premium across delay times
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        'PEP Efficacy Under Source VL Uncertainty — Stochastic Analysis\n'
        'Prevention Theorem | Monte Carlo Extension (n=10,000)',
        fontsize=14, fontweight='bold', y=0.98
    )

    distributions = ['pwid_untreated', 'pwid_mixed_treatment',
                     'general_community', 'acute_infection_enriched']

    # -------------------------------------------------------------------------
    # Panel A: Mean + CI efficacy curves
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    print("Computing stochastic curves...")
    curves = {}
    for dist in distributions:
        model = StochasticPEPModel(vl_distribution=dist, n_simulations=10000)
        curve = model.simulate_timing_curve(max_hours=120, n_timepoints=50)
        curves[dist] = curve
        params = VL_DIST_PARAMS[dist]
        color = params['color']
        hours = curve['hours']

        ax.fill_between(hours,
                        curve['ci_95_lower'] * 100,
                        curve['ci_95_upper'] * 100,
                        alpha=0.12, color=color)
        ax.fill_between(hours,
                        curve['ci_50_lower'] * 100,
                        curve['ci_50_upper'] * 100,
                        alpha=0.25, color=color)
        ax.plot(hours, curve['mean'] * 100,
                color=color, linewidth=2.5, label=params['label'])

    ax.axvline(x=72, color='black', linewidth=1.5, linestyle='--',
               label='72h guideline', alpha=0.7)
    ax.axvline(x=24, color='navy', linewidth=1.5, linestyle='--',
               label='24h optimal', alpha=0.7)
    ax.axhline(y=50, color='gray', linewidth=1, linestyle='-.',
               label='50% threshold', alpha=0.6)

    ax.set_xlabel('Hours from Exposure to PEP Initiation', fontsize=11)
    ax.set_ylabel('Expected PEP Efficacy % (Mean ± 50/95% CI)', fontsize=10)
    ax.set_title('A. Efficacy Under VL Uncertainty\nby Community VL Distribution',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 105)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel B: VL distribution densities
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    vl_range_log = np.linspace(1.3, 7.0, 300)

    for dist in distributions:
        params = VL_DIST_PARAMS[dist]
        density = stats.norm.pdf(vl_range_log,
                                  loc=params['mean_log10'],
                                  scale=params['sd_log10'])
        ax.plot(vl_range_log, density,
                color=params['color'], linewidth=2.5, label=params['label'])
        ax.fill_between(vl_range_log, density, alpha=0.08, color=params['color'])

    for vl_mark, lbl, col in [
        (np.log10(200),    'LLoQ\n200',  'green'),
        (np.log10(10000),  '10K',        'orange'),
        (np.log10(100000), '100K',       'red'),
    ]:
        ax.axvline(x=vl_mark, color=col, linewidth=1, linestyle=':', alpha=0.7)
        ax.text(vl_mark + 0.05, 0.35, lbl, fontsize=8, color=col)

    tick_positions = [1.3, 2, 3, 4, 5, 6, 7]
    tick_labels = ['20\n(LLoD)', '100', '1K', '10K', '100K', '1M', '10M']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlabel('Source Viral Load (copies/mL)', fontsize=11)
    ax.set_ylabel('Probability Density', fontsize=11)
    ax.set_title('B. Community VL Distributions\n(What We Integrate Over)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel C: P(efficacy < 50%) — clinical failure probability
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    for dist in distributions:
        curve = curves[dist]
        params = VL_DIST_PARAMS[dist]
        ax.plot(curve['hours'], curve['p_below_50pct'] * 100,
                color=params['color'], linewidth=2.5, label=params['label'])

    ax.axvline(x=72, color='black', linewidth=1.5, linestyle='--',
               label='72h guideline', alpha=0.7)
    ax.axvline(x=24, color='navy', linewidth=1.5, linestyle='--',
               label='24h optimal', alpha=0.7)
    ax.axhline(y=50, color='gray', linewidth=1, linestyle='-.',
               alpha=0.6, label='50% of patients sub-therapeutic')

    ax.set_xlabel('Hours from Exposure to PEP Initiation', fontsize=11)
    ax.set_ylabel('P(PEP Efficacy < 50%) %', fontsize=11)
    ax.set_title('C. Probability of Sub-Therapeutic PEP\nby Delay and Community',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 105)
    ax.grid(True, color='lightgray', linewidth=0.5)

    # -------------------------------------------------------------------------
    # Panel D: VL knowledge premium
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    delay_times = [0, 3, 6, 9, 12, 18, 24, 36, 48, 60, 72]

    for dist in ['pwid_untreated', 'acute_infection_enriched']:
        params = VL_DIST_PARAMS[dist]
        premium_results = compute_vl_knowledge_premium(
            hours_to_pep_values=delay_times,
            vl_distribution=dist,
            n_simulations=10000
        )
        hours_p = [r['hours'] for r in premium_results]
        premiums = [r['premium_absolute'] * 100 for r in premium_results]

        ax.plot(hours_p, premiums, color=params['color'],
                linewidth=2.5, marker='o', markersize=6,
                label=f"Premium: {params['label'][:30]}")
        ax.fill_between(hours_p, 0, premiums,
                        alpha=0.12, color=params['color'])

    ax.axvline(x=24, color='navy', linewidth=1.5, linestyle='--',
               alpha=0.6, label='24h threshold')
    ax.axhline(y=5, color='gray', linewidth=1, linestyle='-.',
               alpha=0.6, label='5% absolute premium')

    ax.set_xlabel('Hours from Exposure to PEP Initiation', fontsize=11)
    ax.set_ylabel('Efficacy Premium from Knowing Source VL\n(Absolute %, Mean)',
                  fontsize=10)
    ax.set_title('D. Value of Rapid Source VL Testing\n"VL Knowledge Premium"',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 72)
    ax.grid(True, color='lightgray', linewidth=0.5)

    ax.text(0.5, 0.15,
            'Premium = gain in expected efficacy\nif source VL known at PEP initiation',
            transform=ax.transAxes, fontsize=9, ha='center', style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    plt.show()
    return fig, curves


# =============================================================================
# SUMMARY TABLE
# =============================================================================

def print_stochastic_summary():
    """Print clinically interpretable summary table."""
    print("\n" + "=" * 80)
    print("EXPECTED PEP EFFICACY UNDER SOURCE VL UNCERTAINTY")
    print("Mean (95% Credible Interval) | n=10,000 Monte Carlo simulations")
    print("=" * 80)

    timepoints = [0, 12, 24, 48, 72]

    for dist, params in VL_DIST_PARAMS.items():
        model = StochasticPEPModel(vl_distribution=dist, n_simulations=10000)
        print(f"\n{params['label']}")
        mean_row =  "  Mean:    "
        ci_row =    "  95% CI:  "
        futile_row = "  P(futile):"
        for h in timepoints:
            r = model.simulate_at_delay(h)
            mean_row   += f"  {h}h: {r['mean']*100:5.1f}%"
            ci_row     += f"  ({r['ci_95_lower']*100:.0f}-{r['ci_95_upper']*100:.0f}%)"
            futile_row += f"  {r['p_futile']*100:4.1f}%   "
        print(mean_row)
        print(ci_row)
        print(futile_row)

    print("\n" + "=" * 80)
    print("Wide CI = source VL is the dominant unknown.")
    print("Rapid VL testing has highest value when CI is wide AND delay is 12-48h.")
    print("=" * 80)


# =============================================================================
# CSV EXPORT
# =============================================================================

def save_results_csv(output_dir: str = '.'):
    """
    Export all simulation results to CSV files for supplementary data.

    Outputs:
        1. pep_stochastic_efficacy_curves.csv
           Full timing curves (mean + CI bands) for all 4 VL distributions.
           One row per (distribution, timepoint).

        2. pep_stochastic_summary_timepoints.csv
           Summary table at clinically meaningful timepoints (0,12,24,48,72h).
           Includes mean, 95% CI, P(futile) per distribution.

        3. pep_vl_knowledge_premium.csv
           VL knowledge premium at each delay for both PWID distributions.
           Includes raw unknown/known efficacy and absolute + relative premium.

        4. pep_stochastic_ci_width.csv
           95% CI width by timepoint and distribution — captures the
           hourglass shape (regime 1/2/3 boundary identification).
    """
    os.makedirs(output_dir, exist_ok=True)
    distributions = list(VL_DIST_PARAMS.keys())

    # ------------------------------------------------------------------
    # 1. Full timing curves
    # ------------------------------------------------------------------
    print("  Exporting full timing curves...")
    rows_curves = []
    for dist in distributions:
        model = StochasticPEPModel(vl_distribution=dist, n_simulations=10000)
        curve = model.simulate_timing_curve(max_hours=120, n_timepoints=50)
        for i, h in enumerate(curve['hours']):
            rows_curves.append({
                'distribution': dist,
                'distribution_label': VL_DIST_PARAMS[dist]['label'],
                'hours_to_pep': round(h, 2),
                'mean_efficacy': round(curve['mean'][i], 4),
                'median_efficacy': round(curve['median'][i], 4),
                'ci_95_lower': round(curve['ci_95_lower'][i], 4),
                'ci_95_upper': round(curve['ci_95_upper'][i], 4),
                'ci_50_lower': round(curve['ci_50_lower'][i], 4),
                'ci_50_upper': round(curve['ci_50_upper'][i], 4),
                'ci_95_width': round(
                    curve['ci_95_upper'][i] - curve['ci_95_lower'][i], 4),
                'p_below_50pct': round(curve['p_below_50pct'][i], 4),
                'p_below_30pct': round(curve['p_below_30pct'][i], 4),
                'p_futile': round(curve['p_futile'][i], 4),
            })

    df_curves = pd.DataFrame(rows_curves)
    path1 = os.path.join(output_dir, 'pep_stochastic_efficacy_curves.csv')
    df_curves.to_csv(path1, index=False)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
    print(f"    [{timestamp}] Saved: {path1}  ({len(df_curves)} rows)")

    # ------------------------------------------------------------------
    # 2. Summary at key timepoints
    # ------------------------------------------------------------------
    print("  Exporting summary timepoints...")
    timepoints = [0, 6, 12, 18, 24, 36, 48, 72]
    rows_summary = []
    for dist in distributions:
        model = StochasticPEPModel(vl_distribution=dist, n_simulations=10000)
        for h in timepoints:
            r = model.simulate_at_delay(h)
            rows_summary.append({
                'distribution': dist,
                'distribution_label': VL_DIST_PARAMS[dist]['label'],
                'hours_to_pep': h,
                'mean_efficacy_pct': round(r['mean'] * 100, 2),
                'median_efficacy_pct': round(r['median'] * 100, 2),
                'ci_95_lower_pct': round(r['ci_95_lower'] * 100, 2),
                'ci_95_upper_pct': round(r['ci_95_upper'] * 100, 2),
                'ci_50_lower_pct': round(r['ci_50_lower'] * 100, 2),
                'ci_50_upper_pct': round(r['ci_50_upper'] * 100, 2),
                'ci_95_width_pct': round(
                    (r['ci_95_upper'] - r['ci_95_lower']) * 100, 2),
                'p_below_50pct': round(r['p_below_50pct'] * 100, 2),
                'p_below_30pct': round(r['p_below_30pct'] * 100, 2),
                'p_futile_pct': round(r['p_futile'] * 100, 2),
            })

    df_summary = pd.DataFrame(rows_summary)
    path2 = os.path.join(output_dir, 'pep_stochastic_summary_timepoints.csv')
    df_summary.to_csv(path2, index=False)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
    print(f"    [{timestamp}] Saved: {path2}  ({len(df_summary)} rows)")

    # ------------------------------------------------------------------
    # 3. VL knowledge premium — all distributions, fine-grained delays
    # ------------------------------------------------------------------
    print("  Exporting VL knowledge premium...")
    delay_times = [0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 60, 72]
    rows_premium = []
    for dist in distributions:
        results = compute_vl_knowledge_premium(
            hours_to_pep_values=delay_times,
            vl_distribution=dist,
            n_simulations=10000
        )
        for r in results:
            rows_premium.append({
                'distribution': dist,
                'distribution_label': VL_DIST_PARAMS[dist]['label'],
                'hours_to_pep': r['hours'],
                'mean_efficacy_unknown_vl_pct': round(
                    r['mean_efficacy_unknown_vl'] * 100, 2),
                'mean_efficacy_known_vl_pct': round(
                    r['mean_efficacy_known_vl'] * 100, 2),
                'premium_absolute_pct': round(
                    r['premium_absolute'] * 100, 3),
                'premium_relative_pct': round(
                    r['premium_relative'] * 100, 3),
                # Flag non-monotonic regions for methods explanation
                'note': ''
            })

    # Tag the non-monotonic bump in Panel D (18-35h range) —
    # preserving the signal, just labeling it for methods transparency
    df_premium = pd.DataFrame(rows_premium)
    mask_nonmono = (
        (df_premium['hours_to_pep'] >= 18) &
        (df_premium['hours_to_pep'] <= 36) &
        (df_premium['distribution'].isin(
            ['pwid_untreated', 'acute_infection_enriched']))
    )
    df_premium.loc[mask_nonmono, 'note'] = (
        'Non-monotonic region: 2h urgency reduction interacts with '
        'logistic inflection of seeding/integration curves — '
        'preserved as real signal, not smoothed'
    )

    path3 = os.path.join(output_dir, 'pep_vl_knowledge_premium.csv')
    df_premium.to_csv(path3, index=False)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
    print(f"    [{timestamp}] Saved: {path3}  ({len(df_premium)} rows)")

    # ------------------------------------------------------------------
    # 4. CI width by timepoint — regime boundary identification
    # ------------------------------------------------------------------
    print("  Exporting CI width (regime identification)...")
    rows_ciwidth = []
    for dist in distributions:
        # Pull from already-computed curves in df_curves
        subset = df_curves[df_curves['distribution'] == dist].copy()
        # Identify regime boundaries from CI width inflection
        widths = subset['ci_95_width'].values
        hours_arr = subset['hours_to_pep'].values
        peak_idx = int(np.argmax(widths))
        peak_hour = hours_arr[peak_idx]

        for _, row in subset.iterrows():
            h = row['hours_to_pep']
            if h < peak_hour * 0.4:
                regime = '1_early_window'
            elif h <= peak_hour * 1.6:
                regime = '2_critical_window'
            else:
                regime = '3_late_window'
            rows_ciwidth.append({
                'distribution': dist,
                'distribution_label': VL_DIST_PARAMS[dist]['label'],
                'hours_to_pep': row['hours_to_pep'],
                'ci_95_width': row['ci_95_width'],
                'mean_efficacy': row['mean_efficacy'],
                'regime': regime,
                'ci_width_peak_hour': round(peak_hour, 1),
            })

    df_ciwidth = pd.DataFrame(rows_ciwidth)
    path4 = os.path.join(output_dir, 'pep_stochastic_ci_width_regimes.csv')
    df_ciwidth.to_csv(path4, index=False)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
    print(f"    [{timestamp}] Saved: {path4}  ({len(df_ciwidth)} rows)")

    print(f"\n  All CSVs written to: {os.path.abspath(output_dir)}")
    return {
        'efficacy_curves': df_curves,
        'summary_timepoints': df_summary,
        'vl_knowledge_premium': df_premium,
        'ci_width_regimes': df_ciwidth,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')
    print("=" * 80)
    print("PEP STOCHASTIC EXTENSION — VL UNCERTAINTY LAYER")
    print(f"Prevention Theorem | Monte Carlo Analysis | {TIMESTAMP}")
    print("=" * 80)

    print_stochastic_summary()

    print("\nGenerating stochastic figures...")
    fig, curves = plot_stochastic_analysis(
        save_path='pep_stochastic_vl_uncertainty.png'
    )
    print(f"[{TIMESTAMP}] Saved: pep_stochastic_vl_uncertainty.png")

    print("\n" + "=" * 80)
    print("VL KNOWLEDGE PREMIUM — PWID UNTREATED POPULATION")
    print("=" * 80)
    premium_results = compute_vl_knowledge_premium(
        hours_to_pep_values=[0, 6, 12, 24, 48, 72],
        vl_distribution='pwid_untreated'
    )
    for r in premium_results:
        print(f"  {r['hours']:3.0f}h: "
              f"Unknown VL -> {r['mean_efficacy_unknown_vl']*100:.1f}% | "
              f"Known VL -> {r['mean_efficacy_known_vl']*100:.1f}% | "
              f"Premium: +{r['premium_absolute']*100:.1f}pp")

    # --- CSV export ---
    print("\n" + "=" * 80)
    print(f"EXPORTING CSV FILES | {TIMESTAMP}")
    print("=" * 80)
    save_results_csv(output_dir='pep_stochastic_results')

    print("""
THREE REGIMES (from stochastic analysis):

  1. EARLY WINDOW (t < 12h):
     CI narrow regardless of VL distribution.
     Efficacy high even under uncertainty. Route > source VL here.

  2. CRITICAL WINDOW (t = 12-48h):
     CI widens substantially. Source VL is the dominant unknown.
     VL knowledge premium peaks here.
     Rapid POC VL testing has maximum impact in this window.

  3. LATE WINDOW (t > 72h):
     CI narrows again — at LOW efficacy.
     VL irrelevant: PEP is failing regardless.

For PWID under criminalization delays (18-36h), exposures land
squarely in Regime 2. This is the mathematical basis for
point-of-care VL testing as a PWID PEP triage tool.
    """)


if __name__ == "__main__":
    main()