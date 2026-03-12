"""
Combined Prevention Theorem PEP Model: Parenteral, Stochastic, and City-Stratified
================================================================================

This combined script integrates three key components of the Prevention Theorem
PEP modeling framework:

1.  ParenteralExposureModel: Deterministic model for parenteral HIV exposure
    (e.g., PWID, needlestick) with source viral load (VL) as a dynamic window
    compressor. Grounded in Perelson et al. (1996) kinetics.
2.  StochasticPEPModel: Monte Carlo simulation of PEP efficacy under community
    VL uncertainty. Integrates over empirical VL distributions.
3.  City-Stratified Analysis: AIDSVu-grounded city-specific analysis using
    surveillance data (viral suppression, late diagnosis, linkage-to-care) to
    drive city-specific VL distributions and structural delay estimates.

Author: AC Demidont, DO, AAHIVS
Date: March 2026
Updated: 2026-03-10 15:32
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from scipy import stats
from scipy.special import expit
from typing import Dict, List, Tuple, Optional
from datetime import datetime

# Suppress warnings
warnings.filterwarnings('ignore', category=UserWarning)

# --- Path Resolution Utility ---
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

def resolve_project_path(relative_path: str) -> str:
    """Resolve a path relative to the project root."""
    return str(PROJECT_ROOT / relative_path)

# =============================================================================
# CONSTANTS: Perelson 1996 + empirical VL distribution parameters
# =============================================================================

# Log10 VL reference point: median community VL in unsuppressed PWID
# Source: NHBS/NHAS surveillance data; ~30,000 copies/mL in untreated PWID
# NOTE: Perelson 1996 patients had mean plasma VL ~216,000 (log10≈5.3)
REFERENCE_LOG10_VL = 4.5  # ~30,000 copies/mL (NHBS PWID median)

# Seeding midpoint at reference VL (hours post-exposure)
# Parenteral: ~24h vs mucosal ~72h (3x compression from bypassing mucosa)
REFERENCE_SEEDING_MIDPOINT = 24.0

# Integration complete at reference VL (hours post-exposure)
# Grounded in Perelson 1996 Table 2: floor ~22h, outer bound ~62h.
REFERENCE_INTEGRATION_COMPLETE = 56.0

# VL-dependent compression
VL_SEEDING_COMPRESSION_PER_LOG = 4.0
VL_INTEGRATION_COMPRESSION_PER_LOG = 8.0

# Floor: hard biological minimums from Perelson 1996
MIN_SEEDING_MIDPOINT = 6.0       # anchored to virion t1/2 (~6h)
MIN_INTEGRATION_COMPLETE = 22.0  # anchored to eclipse phase (~22h)

# Unsuppressed HIV VL (NHBS PWID median)
VL_UNSUPPRESSED_MEAN_LOG10 = 4.5
VL_UNSUPPRESSED_SD_LOG10 = 1.1

# Suppressed HIV VL on effective ART
VL_SUPPRESSED_MEAN_LOG10 = 1.5
VL_SUPPRESSED_SD_LOG10 = 0.4

# Structural delay mapping: late diagnosis % -> hours delay
LATE_DX_DELAY_SLOPE = 1.2
LATE_DX_DELAY_INTERCEPT = -12.0

# Linkage gap penalty
LINKAGE_GAP_PENALTY_PER_PCT = 0.3

# Inoculum size lookup by exposure type
INOCULUM_BY_EXPOSURE = {
    'needlestick_hollow':    150,   # Hollow-bore needle, ~1mL blood
    'needlestick_solid':     10,    # Solid suture needle
    'pwid_shared_needle':    500,   # PWID: residual blood in syringe
    'pwid_shared_cooker':    50,    # Shared cooker/water
    'blood_transfusion':     1e8,   # Transfusion: massive inoculum
}

# VL Distribution Scenario Parameters
VL_DIST_PARAMS = {
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
    'acute_infection_enriched': {
        'mean_log10': 5.0,
        'sd_log10': 0.9,
        'label': 'Acute-infection-enriched network',
        'color': 'darkred'
    }
}

# AIDSVu IDU transmission columns
COL_MALE_IDU_PREV = '2023 Male\nInjection\nDrug Use\nPrevalence\nPercent'
COL_FEMALE_IDU_PREV = '2023\nFemale\nInjection\nDrug Use\nPrevalence\nPercent'
COL_MALE_IDU_DX = '2023 Male\nInjection\nDrug Use\nNew\nDiagnoses\nPercent'
COL_FEMALE_IDU_DX = '2023\nFemale\nInjection\nDrug Use\nNew\nDiagnoses\nPercent'
COL_VS_PCT = '2023 Viral\nSuppression\nPercent'
COL_LATE_DX_PCT = '2023 Late\nDiagnoses\nPercent'
COL_LINKAGE_PCT = '2023\nLinkage\nto Care\nPercent'
COL_NEW_DX_RATE = '2023 New\nDiagnoses\nRate'
COL_NEW_DX_CASES = '2023 New\nDiagnoses\nCases'
COL_PREVALENCE_RATE = '2023\nPrevalence\nRate'
COL_PREVALENCE_CASES = '2023\nPrevalence\nCases'
COL_GINI = '2023 Gini\nCoefficient'
COL_POVERTY_PCT = '2023\nPercent\nLiving\nin\nPoverty'
COL_HOUSING_BURDEN = '2023\nPercent\nLiving\nwith\nSevere\nHousing\nCost\nBurden'
COL_REGION = 'Region'
COL_STATE = 'State'
COL_POPULATION = '2020 City\nCensus\nPopulation'

# =============================================================================
# DETERMINISTIC MODEL
# =============================================================================

class ParenteralExposureModel:
    """
    Models PEP efficacy for parenteral HIV exposure with source VL integration.
    """
    def __init__(self,
                 source_viral_load: float = 30000.0,
                 exposure_type: str = 'pwid_shared_needle'):
        self.source_viral_load = max(source_viral_load, 1.0)
        self.exposure_type = exposure_type
        self.log10_vl = np.log10(self.source_viral_load)

        # Parenteral timing (no mucosal phase)
        self.pep_onset_hours = 2.0
        self.pep_efficacy_peak = 0.995
        self.seeding_midpoint = self._compute_seeding_midpoint()
        self.integration_complete = self._compute_integration_complete()
        self.inoculum_virions = INOCULUM_BY_EXPOSURE.get(exposure_type, 150)

    def _compute_seeding_midpoint(self) -> float:
        delta_log = self.log10_vl - REFERENCE_LOG10_VL
        compressed = REFERENCE_SEEDING_MIDPOINT - (delta_log * VL_SEEDING_COMPRESSION_PER_LOG)
        return max(compressed, MIN_SEEDING_MIDPOINT)

    def _compute_integration_complete(self) -> float:
        delta_log = self.log10_vl - REFERENCE_LOG10_VL
        compressed = REFERENCE_INTEGRATION_COMPLETE - (delta_log * VL_INTEGRATION_COMPRESSION_PER_LOG)
        return max(compressed, MIN_INTEGRATION_COMPLETE)

    def effective_pep_window(self) -> float:
        return max(self.integration_complete - self.pep_onset_hours, 0.0)

    def p_seeding_initiated(self, hours_post_exposure: float) -> float:
        k = 0.15
        return float(expit(k * (hours_post_exposure - self.seeding_midpoint)))

    def p_integration_complete(self, hours_post_exposure: float) -> float:
        k = 0.20
        return float(expit(k * (hours_post_exposure - self.integration_complete)))

    def pep_efficacy(self, hours_to_pep: float, adherence: float = 1.0) -> Dict:
        effective_time = hours_to_pep + self.pep_onset_hours
        p_seeded = self.p_seeding_initiated(effective_time)
        p_integrated = self.p_integration_complete(effective_time)

        eff_not_seeded = self.pep_efficacy_peak * adherence
        eff_seeded_not_integrated = 0.40 * adherence
        eff_integrated = 0.0

        overall_efficacy = (
            (1 - p_seeded) * eff_not_seeded +
            (p_seeded - p_integrated) * eff_seeded_not_integrated +
            p_integrated * eff_integrated
        )

        return {
            'hours_to_pep': hours_to_pep,
            'source_viral_load': self.source_viral_load,
            'seeding_midpoint': self.seeding_midpoint,
            'integration_complete': self.integration_complete,
            'pep_efficacy': overall_efficacy,
            'p_seeding_initiated': p_seeded,
            'p_integration_complete': p_integrated,
        }

    def simulate_timing_curve(self, max_hours: float = 120.0, n_points: int = 200) -> Dict:
        hours = np.linspace(0, max_hours, n_points)
        effs = [self.pep_efficacy(h)['pep_efficacy'] for h in hours]
        return {'hours': hours, 'efficacy': np.array(effs)}

# =============================================================================
# STOCHASTIC MODEL
# =============================================================================

class StochasticPEPModel:
    """
    Monte Carlo simulation of PEP efficacy under source VL uncertainty.
    """
    def __init__(self,
                 vl_distribution: str = 'pwid_untreated',
                 exposure_type: str = 'pwid_shared_needle',
                 n_simulations: int = 10000,
                 random_seed: int = 42):
        self.vl_dist_params = VL_DIST_PARAMS[vl_distribution]
        self.exposure_type = exposure_type
        self.n_simulations = n_simulations
        self.rng = np.random.default_rng(random_seed)

        # Pre-sample VL values
        sampled_log10 = self.rng.normal(loc=self.vl_dist_params['mean_log10'],
                                         scale=self.vl_dist_params['sd_log10'],
                                         size=n_simulations)
        self.sampled_log10_vl = np.clip(sampled_log10, 1.3, 7.0)
        self.sampled_vl = 10 ** self.sampled_log10_vl

    def simulate_at_delay(self, hours_to_pep: float) -> Dict:
        """Run Monte Carlo at a single delay time."""
        efficacies = []
        for vl in self.sampled_vl:
            model = ParenteralExposureModel(source_viral_load=vl, exposure_type=self.exposure_type)
            efficacies.append(model.pep_efficacy(hours_to_pep)['pep_efficacy'])
        
        efficacies = np.array(efficacies)
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
            'p_futile': np.mean(efficacies < 0.1)
        }

    def simulate_timing_curve(self, max_hours: float = 120.0, n_timepoints: int = 50) -> Dict:
        """Run Monte Carlo across full timing curve."""
        hours = np.linspace(0, max_hours, n_timepoints)
        results = {k: [] for k in ['hours', 'mean', 'median', 'ci_95_lower', 'ci_95_upper', 
                                   'ci_50_lower', 'ci_50_upper', 'p_below_50pct', 'p_futile']}
        results['hours'] = hours
        for h in hours:
            r = self.simulate_at_delay(h)
            for key in results:
                if key != 'hours': results[key].append(r[key])
        for key in results:
            if key != 'hours': results[key] = np.array(results[key])
        return results

# =============================================================================
# VISUALIZATION & UTILITIES
# =============================================================================

def compute_vl_knowledge_premium(hours: List[float], vl_dist: str = 'pwid_untreated') -> List[Dict]:
    params = VL_DIST_PARAMS[vl_dist]
    rng = np.random.default_rng(42)
    sampled_vl = 10 ** np.clip(rng.normal(params['mean_log10'], params['sd_log10'], 10000), 1.3, 7.0)
    
    results = []
    for h in hours:
        effs_unk, effs_knw = [], []
        for vl in sampled_vl:
            effs_unk.append(ParenteralExposureModel(source_viral_load=vl).pep_efficacy(h)['pep_efficacy'])
            effs_knw.append(ParenteralExposureModel(source_viral_load=vl).pep_efficacy(max(h-2.0, 0))['pep_efficacy'] if vl > 10000 else effs_unk[-1])
        results.append({'hours': h, 'premium': np.mean(effs_knw) - np.mean(effs_unk)})
    return results

def plot_stochastic_analysis(save_path: str):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    dists = list(VL_DIST_PARAMS.keys())
    
    # A. Efficacy Curves
    ax = axes[0, 0]
    for d in dists:
        curve = StochasticPEPModel(vl_distribution=d).simulate_timing_curve()
        ax.plot(curve['hours'], curve['mean']*100, label=VL_DIST_PARAMS[d]['label'], color=VL_DIST_PARAMS[d]['color'])
        ax.fill_between(curve['hours'], curve['ci_95_lower']*100, curve['ci_95_upper']*100, alpha=0.1, color=VL_DIST_PARAMS[d]['color'])
    ax.set_title("A. Efficacy Under VL Uncertainty"); ax.legend(); ax.set_xlabel("Hours"); ax.set_ylabel("Efficacy %")

    # B. VL Distributions
    ax = axes[0, 1]
    x = np.linspace(1.3, 7, 200)
    for d in dists:
        p = VL_DIST_PARAMS[d]
        ax.plot(x, stats.norm.pdf(x, p['mean_log10'], p['sd_log10']), label=p['label'], color=p['color'])
    ax.set_title("B. Source VL Distributions"); ax.legend()

    # C. Failure Prob
    ax = axes[1, 0]
    for d in dists:
        curve = StochasticPEPModel(vl_distribution=d).simulate_timing_curve()
        ax.plot(curve['hours'], curve['p_below_50pct']*100, label=VL_DIST_PARAMS[d]['label'], color=VL_DIST_PARAMS[d]['color'])
    ax.set_title("C. P(Efficacy < 50%)"); ax.legend()

    # D. Knowledge Premium
    ax = axes[1, 1]
    for d in ['pwid_untreated', 'acute_infection_enriched']:
        prem = compute_vl_knowledge_premium([0, 6, 12, 18, 24, 36, 48, 72], vl_dist=d)
        ax.plot([r['hours'] for r in prem], [r['premium']*100 for r in prem], label=VL_DIST_PARAMS[d]['label'], marker='o', color=VL_DIST_PARAMS[d]['color'])
    ax.set_title("D. VL Knowledge Premium"); ax.legend()

    plt.tight_layout()
    plt.savefig(save_path); plt.close()

def plot_vl_sweep(save_path: str):
    fig, ax = plt.subplots(figsize=(10, 6))
    vls = [200, 1000, 10000, 100000, 500000]
    cmap = plt.cm.RdYlGn_r
    for vl in vls:
        curve = ParenteralExposureModel(source_viral_load=vl).simulate_timing_curve()
        ax.plot(curve['hours'], curve['efficacy']*100, label=f"VL {vl:,}", color=cmap(np.log10(vl)/7))
    ax.set_title("Parenteral PEP Efficacy by Source VL"); ax.legend(); ax.set_xlabel("Hours"); ax.set_ylabel("Efficacy %")
    plt.savefig(save_path); plt.close()

# =============================================================================
# CITY PROFILE PARSER & DERIVATION
# =============================================================================

def parse_city_profile(filepath: str) -> Optional[Dict]:
    city_name = os.path.basename(filepath).split('_AIDSVu')[0]
    try:
        df = pd.read_excel(filepath, sheet_name='City Profile', header=None)
        headers = list(df.iloc[3].values)
        data = list(df.iloc[4].values)
        lookup = {str(h): v for h, v in zip(headers, data) if pd.notna(h)}

        def get(col, default=np.nan):
            v = lookup.get(col, default)
            if v == -1 or v == '-1': return np.nan
            try: return float(v)
            except: return default

        male_idu_prev = get(COL_MALE_IDU_PREV)
        female_idu_prev = get(COL_FEMALE_IDU_PREV)
        male_pct = get('2023 Male\nPrevalence\nPercent', 75) / 100
        female_pct = 1 - male_pct
        idu_prev = male_idu_prev * male_pct + female_idu_prev * female_pct if pd.notna(male_idu_prev) else np.nan

        return {
            'city': city_name,
            'state': str(lookup.get(COL_STATE, '')),
            'region': str(lookup.get(COL_REGION, '')),
            'viral_suppression_pct': get(COL_VS_PCT),
            'late_dx_pct': get(COL_LATE_DX_PCT),
            'linkage_to_care_pct': get(COL_LINKAGE_PCT),
            'idu_prevalence_pct': idu_prev,
            'idu_new_dx_pct': (get(COL_MALE_IDU_DX)*male_pct + get(COL_FEMALE_IDU_DX)*female_pct) if pd.notna(get(COL_MALE_IDU_DX)) else np.nan,
            'new_dx_cases': get(COL_NEW_DX_CASES),
            'prevalence_rate': get(COL_PREVALENCE_RATE),
            'gini': get(COL_GINI),
            'poverty_pct': get(COL_POVERTY_PCT),
        }
    except Exception as e:
        print(f"  WARNING: Failed to parse {city_name}: {e}")
        return None

def derive_vl_distribution(vs_pct: float, idu_prev: float = np.nan) -> Dict:
    vs_frac = vs_pct / 100.0
    if pd.notna(idu_prev) and idu_prev > 0:
        vs_frac = max(vs_frac - min(idu_prev * 0.0012, 0.12), 0.1)
    
    unsupp_frac = 1.0 - vs_frac
    mean_log10 = vs_frac * VL_SUPPRESSED_MEAN_LOG10 + unsupp_frac * VL_UNSUPPRESSED_MEAN_LOG10
    
    e_x2_supp = VL_SUPPRESSED_SD_LOG10**2 + VL_SUPPRESSED_MEAN_LOG10**2
    e_x2_unsupp = VL_UNSUPPRESSED_SD_LOG10**2 + VL_UNSUPPRESSED_MEAN_LOG10**2
    var_mix = (vs_frac * e_x2_supp + unsupp_frac * e_x2_unsupp) - mean_log10**2
    sd_log10 = max(np.sqrt(var_mix), 0.3)

    return {'mean_log10': round(mean_log10, 3), 'sd_log10': round(sd_log10, 3)}

def derive_structural_delay(late_dx: float, linkage: float, gini: float, poverty: float) -> Dict:
    late_dx_h = max(LATE_DX_DELAY_SLOPE * (late_dx - 10.0) + LATE_DX_DELAY_INTERCEPT, 0.0)
    linkage_h = max(90.0 - linkage, 0.0) * LINKAGE_GAP_PENALTY_PER_PCT
    sdoh_h = 0.0
    if pd.notna(gini) and pd.notna(poverty):
        sdoh_h = min(max(gini - 0.40, 0) * 20.0 + max(poverty - 15.0, 0) * 0.1, 6.0)
    
    total = late_dx_h + linkage_h + sdoh_h
    cat = 'low_barrier' if total < 8 else ('moderate_barrier' if total < 18 else ('high_barrier' if total < 30 else 'severe_barrier'))
    return {'estimated_structural_delay_h': round(total, 1), 'delay_category': cat}

def build_city_profiles(data_dir: str = None) -> pd.DataFrame:
    if data_dir is None: data_dir = resolve_project_path('aidsvu datasets')
    files = sorted(glob.glob(os.path.join(data_dir, '*AIDSVu*.xlsx')))
    if not files: print("ERROR: No AIDSVu files found."); return pd.DataFrame()
    
    rows = []
    for fp in files:
        profile = parse_city_profile(fp)
        if profile:
            vl = derive_vl_distribution(profile['viral_suppression_pct'], profile['idu_prevalence_pct'])
            delay = derive_structural_delay(profile['late_dx_pct'], profile['linkage_to_care_pct'], profile['gini'], profile['poverty_pct'])
            rows.append({**profile, **vl, **delay})
    return pd.DataFrame(rows)

# =============================================================================
# ANALYSIS PIPELINE
# =============================================================================

def run_stratified_analysis(city_df: pd.DataFrame, n_sims: int = 5000) -> pd.DataFrame:
    rng = np.random.default_rng(42)
    results = []
    for _, city in city_df.iterrows():
        delay_h = city['estimated_structural_delay_h']
        sampled_log10 = np.clip(rng.normal(city['mean_log10'], city['sd_log10'], n_sims), 1.3, 7.0)
        sampled_vl = 10**sampled_log10
        
        effs = []
        for vl in sampled_vl:
            effs.append(ParenteralExposureModel(source_viral_load=vl).pep_efficacy(delay_h)['pep_efficacy'])
        
        effs = np.array(effs)
        results.append({
            'city': city['city'],
            'structural_delay_h': delay_h,
            'pep_mean_efficacy_pct': round(np.mean(effs)*100, 2),
            'p_futile_pct': round(np.mean(effs < 0.1)*100, 2)
        })
    return pd.DataFrame(results).sort_values('pep_mean_efficacy_pct')

# =============================================================================
# MAIN
# =============================================================================

def main():
    TIMESTAMP = datetime.now().strftime('%Y-%m-%d %H:%M')
    print("=" * 80)
    print("PREVENTION THEOREM: COMBINED STOCHASTIC PARENTERAL CITY-STRATIFIED MODEL")
    print(f"[{TIMESTAMP}]")
    print("=" * 80)

    output_dir = resolve_project_path('results/city_analysis')
    os.makedirs(output_dir, exist_ok=True)

    # 1. City Profiles
    city_df = build_city_profiles()
    if city_df.empty: return
    city_df.to_csv(os.path.join(output_dir, 'city_profiles_combined.csv'), index=False)
    print(f"[{TIMESTAMP}] Saved: city_profiles_combined.csv ({len(city_df)} cities)")

    # 2. Stratified Analysis
    print("Running stratified stochastic analysis...")
    results = run_stratified_analysis(city_df)
    results.to_csv(os.path.join(output_dir, 'city_pep_efficacy_combined.csv'), index=False)
    print(f"[{TIMESTAMP}] Saved: city_pep_efficacy_combined.csv")
    
    print("\nTop 10 Cities by PEP Efficacy (Ranked worst to best):")
    print(results.head(10).to_string(index=False))

    # 3. Visualizations
    print("Generating figures...")
    path_fig1 = os.path.join(output_dir, 'pep_stochastic_vl_uncertainty_combined.png')
    plot_stochastic_analysis(path_fig1)
    print(f"[{TIMESTAMP}] Saved: {path_fig1}")

    path_fig2 = os.path.join(output_dir, 'pep_parenteral_vl_sweep_combined.png')
    plot_vl_sweep(path_fig2)
    print(f"[{TIMESTAMP}] Saved: {path_fig2}")

    # 4. CSV Exports
    print("Exporting results to CSV...")
    # Summary Table (key timepoints)
    timepoints = [0, 12, 24, 48, 72]
    summary_rows = []
    for d_name in VL_DIST_PARAMS:
        model = StochasticPEPModel(vl_distribution=d_name)
        for tp in timepoints:
            r = model.simulate_at_delay(tp)
            summary_rows.append({
                'distribution': d_name,
                'hours_to_pep': tp,
                'mean_efficacy_pct': round(r['mean']*100, 2),
                'ci_95_lower_pct': round(r['ci_95_lower']*100, 2),
                'ci_95_upper_pct': round(r['ci_95_upper']*100, 2),
                'p_futile_pct': round(r['p_futile']*100, 2)
            })
    path_summary = os.path.join(output_dir, 'pep_stochastic_summary_combined.csv')
    pd.DataFrame(summary_rows).to_csv(path_summary, index=False)
    print(f"[{TIMESTAMP}] Saved: {path_summary}")

    print(f"\n[{TIMESTAMP}] Analysis complete. Outputs saved to results/city_analysis/")

if __name__ == "__main__":
    main()
