"""
Regenerate all PEP analysis CSVs using Perelson-corrected model constants.

Changes from previous run:
  - REFERENCE_INTEGRATION_COMPLETE: 60h → 56h
  - MIN_SEEDING_MIDPOINT:           4h  → 6h  (Perelson virion t½ ~6h)
  - MIN_INTEGRATION_COMPLETE:       12h → 22h (Perelson eclipse phase ~22h)
  - Sign convention fix: efficacy_lost column renamed and sign corrected
  - VL citation fix: log10=4.5 now correctly attributed to NHBS, not Perelson

Perelson AS et al. Science 1996;271(5255):1582. PMID: 8599114
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy import stats
import os

# =============================================================================
# CORRECTED CONSTANTS (Perelson-grounded)
# =============================================================================

REFERENCE_LOG10_VL              = 4.5   # NHBS PWID median
REFERENCE_SEEDING_MIDPOINT      = 24.0  # hours
REFERENCE_INTEGRATION_COMPLETE  = 56.0  # hours (corrected from 60)
VL_SEEDING_COMPRESSION_PER_LOG  = 4.0
VL_INTEGRATION_COMPRESSION_PER_LOG = 8.0
MIN_SEEDING_MIDPOINT             = 6.0  # hours — Perelson virion t½ = 0.24d = 5.8h
MIN_INTEGRATION_COMPLETE         = 22.0 # hours — Perelson eclipse phase = 0.9d = 21.6h
PEP_ONSET_HOURS                  = 2.0
PEP_EFFICACY_PEAK                = 0.995


# =============================================================================
# CORE MODEL FUNCTIONS
# =============================================================================

def compute_integration_complete(log10_vl):
    delta = log10_vl - REFERENCE_LOG10_VL
    val = REFERENCE_INTEGRATION_COMPLETE - (delta * VL_INTEGRATION_COMPRESSION_PER_LOG)
    return max(val, MIN_INTEGRATION_COMPLETE)

def compute_seeding_midpoint(log10_vl):
    delta = log10_vl - REFERENCE_LOG10_VL
    val = REFERENCE_SEEDING_MIDPOINT - (delta * VL_SEEDING_COMPRESSION_PER_LOG)
    return max(val, MIN_SEEDING_MIDPOINT)

def pep_efficacy_single(hours_to_pep, log10_vl):
    seeding_mid = compute_seeding_midpoint(log10_vl)
    int_complete = compute_integration_complete(log10_vl)
    eff_t = hours_to_pep + PEP_ONSET_HOURS
    k_s, k_i = 0.15, 0.20
    p_s = float(expit(k_s * (eff_t - seeding_mid)))
    p_i = float(expit(k_i * (eff_t - int_complete)))
    return (1 - p_s) * PEP_EFFICACY_PEAK + (p_s - p_i) * 0.5 + p_i * 0.0

def simulate_stochastic(mean_log10, sd_log10, hours, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    vl_samples = np.clip(rng.normal(mean_log10, sd_log10, n), 1.3, 7.0)
    effs = np.array([pep_efficacy_single(hours, vl) for vl in vl_samples])
    return {
        'mean': np.mean(effs),
        'median': np.median(effs),
        'ci95_lo': np.percentile(effs, 2.5),
        'ci95_hi': np.percentile(effs, 97.5),
        'ci50_lo': np.percentile(effs, 25),
        'ci50_hi': np.percentile(effs, 75),
        'ci95_width': np.percentile(effs, 97.5) - np.percentile(effs, 2.5),
        'p_futile': np.mean(effs < 0.1),
        'p_below_50': np.mean(effs < 0.5),
        'p_below_30': np.mean(effs < 0.3),
    }


# =============================================================================
# VL DISTRIBUTIONS
# =============================================================================

VL_DIST_PARAMS = {
    'pwid_untreated':            {'mean_log10': 4.5, 'sd_log10': 1.1,
                                  'label': 'PWID (untreated/unsuppressed)'},
    'pwid_mixed_treatment':      {'mean_log10': 3.2, 'sd_log10': 1.5,
                                  'label': 'PWID (mixed treatment status)'},
    'general_community':         {'mean_log10': 3.8, 'sd_log10': 1.2,
                                  'label': 'General HIV+ community'},
    'acute_infection_enriched':  {'mean_log10': 5.0, 'sd_log10': 0.9,
                                  'label': 'Acute-infection-enriched network'},
}

TIMEPOINTS = [0, 6, 12, 18, 24, 36, 48, 72]
TIMEPOINTS_DENSE = list(np.linspace(0, 120, 50))


# =============================================================================
# 1. STOCHASTIC SUMMARY TIMEPOINTS CSV
# =============================================================================

print("Generating stochastic_summary_timepoints...")
rows = []
for dist, params in VL_DIST_PARAMS.items():
    for h in TIMEPOINTS:
        r = simulate_stochastic(params['mean_log10'], params['sd_log10'], h)
        rows.append({
            'distribution': dist,
            'distribution_label': params['label'],
            'hours_to_pep': h,
            'mean_efficacy_pct': round(r['mean']*100, 2),
            'median_efficacy_pct': round(r['median']*100, 2),
            'ci_95_lower_pct': round(r['ci95_lo']*100, 2),
            'ci_95_upper_pct': round(r['ci95_hi']*100, 2),
            'ci_50_lower_pct': round(r['ci50_lo']*100, 2),
            'ci_50_upper_pct': round(r['ci50_hi']*100, 2),
            'ci_95_width_pct': round(r['ci95_width']*100, 2),
            'p_below_50pct': round(r['p_below_50']*100, 2),
            'p_below_30pct': round(r['p_below_30']*100, 2),
            'p_futile_pct': round(r['p_futile']*100, 2),
        })
df_summary = pd.DataFrame(rows)
df_summary.to_csv('/mnt/user-data/outputs/pep_stochastic_summary_timepoints_corrected.csv', index=False)
print(f"  Done. {len(df_summary)} rows.")


# =============================================================================
# 2. EFFICACY CURVES CSV (50 timepoints)
# =============================================================================

print("Generating efficacy_curves...")
rows = []
for dist, params in VL_DIST_PARAMS.items():
    for h in TIMEPOINTS_DENSE:
        r = simulate_stochastic(params['mean_log10'], params['sd_log10'], h, n=5000)
        rows.append({
            'distribution': dist,
            'distribution_label': params['label'],
            'hours_to_pep': round(h, 2),
            'mean_efficacy': round(r['mean'], 4),
            'median_efficacy': round(r['median'], 4),
            'ci_95_lower': round(r['ci95_lo'], 4),
            'ci_95_upper': round(r['ci95_hi'], 4),
            'ci_50_lower': round(r['ci50_lo'], 4),
            'ci_50_upper': round(r['ci50_hi'], 4),
            'ci_95_width': round(r['ci95_width'], 4),
            'p_below_50pct': round(r['p_below_50'], 4),
            'p_below_30pct': round(r['p_below_30'], 4),
            'p_futile': round(r['p_futile'], 4),
        })
df_curves = pd.DataFrame(rows)
df_curves.to_csv('/mnt/user-data/outputs/pep_stochastic_efficacy_curves_corrected.csv', index=False)
print(f"  Done. {len(df_curves)} rows.")


# =============================================================================
# 3. CI WIDTH REGIMES CSV
# =============================================================================

print("Generating ci_width_regimes...")
rows = []
for dist, params in VL_DIST_PARAMS.items():
    peak_ci_h = None
    peak_ci_w = 0
    # Find peak CI width first
    for h in TIMEPOINTS_DENSE:
        r = simulate_stochastic(params['mean_log10'], params['sd_log10'], h, n=5000)
        if r['ci95_width'] > peak_ci_w:
            peak_ci_w = r['ci95_width']
            peak_ci_h = h

    for h in TIMEPOINTS_DENSE:
        r = simulate_stochastic(params['mean_log10'], params['sd_log10'], h, n=5000)
        if h < 22:
            regime = '1_early_window'
        elif h <= 88:
            regime = '2_critical_window'
        else:
            regime = '3_late_window'
        rows.append({
            'distribution': dist,
            'distribution_label': params['label'],
            'hours_to_pep': round(h, 2),
            'ci_95_width': round(r['ci95_width'], 4),
            'mean_efficacy': round(r['mean'], 4),
            'regime': regime,
            'ci_width_peak_hour': round(peak_ci_h, 1),
        })

df_regimes = pd.DataFrame(rows)
df_regimes.to_csv('/mnt/user-data/outputs/pep_stochastic_ci_width_regimes_corrected.csv', index=False)
print(f"  Done. {len(df_regimes)} rows.")


# =============================================================================
# 4. VL KNOWLEDGE PREMIUM CSV
# =============================================================================

print("Generating vl_knowledge_premium...")
PREMIUM_TIMEPOINTS = [0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 60, 72]
NON_MONO_RANGE = (18, 36)

rows = []
for dist, params in VL_DIST_PARAMS.items():
    rng = np.random.default_rng(42)
    vl_samples = np.clip(rng.normal(params['mean_log10'], params['sd_log10'], 10000), 1.3, 7.0)

    for h in PREMIUM_TIMEPOINTS:
        # Unknown VL: mean over full distribution
        effs_unknown = np.array([pep_efficacy_single(h, vl) for vl in vl_samples])
        mean_unknown = np.mean(effs_unknown)

        # Known VL: high-VL cases get 2h urgency reduction
        effs_known = []
        for vl in vl_samples:
            if vl > np.log10(10000):  # high VL known → 2h urgency reduction
                eff = pep_efficacy_single(max(h - 2, 0), vl)
            else:
                eff = pep_efficacy_single(h, vl)
            effs_known.append(eff)
        mean_known = np.mean(effs_known)
        premium_abs = mean_known - mean_unknown
        premium_rel = (premium_abs / mean_unknown * 100) if mean_unknown > 0 else 0

        note = ''
        if NON_MONO_RANGE[0] <= h <= NON_MONO_RANGE[1]:
            note = ('Non-monotonic region: 2h urgency reduction interacts with '
                    'logistic inflection of seeding/integration curves — '
                    'preserved as real signal, not smoothed')

        rows.append({
            'distribution': dist,
            'distribution_label': params['label'],
            'hours_to_pep': h,
            'mean_efficacy_unknown_vl_pct': round(mean_unknown * 100, 2),
            'mean_efficacy_known_vl_pct': round(mean_known * 100, 2),
            'premium_absolute_pct': round(premium_abs * 100, 3),
            'premium_relative_pct': round(premium_rel, 3),
            'note': note,
        })

df_premium = pd.DataFrame(rows)
df_premium.to_csv('/mnt/user-data/outputs/pep_vl_knowledge_premium_corrected.csv', index=False)
print(f"  Done. {len(df_premium)} rows.")


# =============================================================================
# 5. CITY EFFICACY CSVs — Regenerate from corrected model
# =============================================================================

print("Regenerating city efficacy CSVs...")
city_vl = pd.read_csv('/mnt/user-data/uploads/1773169285545_city_vl_profiles.csv')

city_rows = []
for _, row in city_vl.iterrows():
    delay = row['estimated_structural_delay_h']
    mean_log = row['mean_log10']
    sd_log   = row['sd_log10']

    r = simulate_stochastic(mean_log, sd_log, delay, n=10000)

    city_rows.append({
        'city': row['city'],
        'state': row['state'],
        'region': row['region'],
        'viral_suppression_pct': row['viral_suppression_pct'],
        'late_dx_pct': row['late_dx_pct'],
        'linkage_to_care_pct': row['linkage_to_care_pct'],
        'idu_prevalence_pct': row['idu_prevalence_pct'],
        'vl_mean_log10': mean_log,
        'vl_sd_log10': sd_log,
        'vl_profile': row['vl_profile'],
        'structural_delay_h': delay,
        'delay_category': row['delay_category'],
        'pep_mean_efficacy_pct': round(r['mean']*100, 2),
        'pep_median_efficacy_pct': round(r['median']*100, 2),
        'pep_ci95_lower_pct': round(r['ci95_lo']*100, 2),
        'pep_ci95_upper_pct': round(r['ci95_hi']*100, 2),
        'p_futile_pct': round(r['p_futile']*100, 2),
        'p_below_50pct': round(r['p_below_50']*100, 2),
        'estimated_pwid_exposures': row['estimated_annual_pwid_exposures'],
        'new_dx_cases': row['new_dx_cases'],
        'prevalence_rate': row['prevalence_rate'],
        'gini': row['gini'],
        'poverty_pct': row['poverty_pct'],
    })

df_city = pd.DataFrame(city_rows).sort_values('pep_mean_efficacy_pct')
df_city.to_csv('/mnt/user-data/outputs/city_pep_efficacy_results_corrected.csv', index=False)
print(f"  City efficacy: {len(df_city)} cities.")

# Counterfactual at 24h — fixed sign convention
cf_rows = []
for _, row in city_vl.iterrows():
    mean_log = row['mean_log10']
    sd_log   = row['sd_log10']
    r_24h = simulate_stochastic(mean_log, sd_log, 24.0, n=10000)
    cf_rows.append({
        'city': row['city'],
        'state': row['state'],
        'region': row['region'],
        'vl_mean_log10': mean_log,
        'structural_delay_h': row['estimated_structural_delay_h'],
        'pep_mean_efficacy_pct_24h_counterfactual': round(r_24h['mean']*100, 2),
        'pep_ci95_lower_pct': round(r_24h['ci95_lo']*100, 2),
        'pep_ci95_upper_pct': round(r_24h['ci95_hi']*100, 2),
    })

df_cf = pd.DataFrame(cf_rows)
df_city_cf = df_city.merge(df_cf[['city','pep_mean_efficacy_pct_24h_counterfactual',
                                   'pep_ci95_lower_pct','pep_ci95_upper_pct']]
                            .rename(columns={'pep_ci95_lower_pct':'cf_ci95_lower_pct',
                                            'pep_ci95_upper_pct':'cf_ci95_upper_pct'}),
                           on='city')
df_city_cf.to_csv('/mnt/user-data/outputs/city_pep_efficacy_24h_counterfactual_corrected.csv', index=False)
print(f"  Counterfactual: {len(df_city_cf)} cities.")

# Structural delay cost — FIXED sign convention
# efficacy_gain_vs_24h_counterfactual_pp = actual - counterfactual
# Positive = city benefits from lower structural delay (vs 24h baseline)
# Negative = city WORSE than 24h (only Hartford expected)
cost_rows = []
for _, row in df_city_cf.iterrows():
    gain = row['pep_mean_efficacy_pct'] - row['pep_mean_efficacy_pct_24h_counterfactual']
    cost_rows.append({
        'city': row['city'],
        'structural_delay_h': row['structural_delay_h'],
        'pep_mean_efficacy_pct': row['pep_mean_efficacy_pct'],
        'pep_24h_counterfactual_pct': row['pep_mean_efficacy_pct_24h_counterfactual'],
        'efficacy_gain_vs_24h_counterfactual_pp': round(gain, 2),
        # Note: positive = better than 24h; negative = worse than 24h
    })

df_cost = pd.DataFrame(cost_rows).sort_values('efficacy_gain_vs_24h_counterfactual_pp')
df_cost.to_csv('/mnt/user-data/outputs/city_structural_delay_cost_corrected.csv', index=False)
print(f"  Delay cost: {len(df_cost)} cities.")

# =============================================================================
# SUMMARY OF CHANGES
# =============================================================================

print("\n" + "="*65)
print("SUMMARY: KEY NUMERICAL CHANGES FROM CORRECTION")
print("="*65)

# PWID untreated at 24h before and after
old_24h = 65.46  # from original CSV
new_24h = df_summary[
    (df_summary['distribution']=='pwid_untreated') &
    (df_summary['hours_to_pep']==24)
]['mean_efficacy_pct'].values[0]
print(f"\nPWID untreated at 24h:")
print(f"  Original: {old_24h:.1f}%")
print(f"  Corrected: {new_24h:.1f}%")
print(f"  Delta: {new_24h-old_24h:+.1f}pp")

# Hartford
old_hartford = 78.84
new_hartford = df_city[df_city['city']=='Hartford']['pep_mean_efficacy_pct'].values[0]
print(f"\nHartford (24.4h structural delay):")
print(f"  Original: {old_hartford:.2f}%")
print(f"  Corrected: {new_hartford:.2f}%")
print(f"  Delta: {new_hartford-old_hartford:+.2f}pp")

# Milwaukee
old_milw = 98.34
new_milw = df_city[df_city['city']=='Milwaukee']['pep_mean_efficacy_pct'].values[0]
print(f"\nMilwaukee (1.9h structural delay):")
print(f"  Original: {old_milw:.2f}%")
print(f"  Corrected: {new_milw:.2f}%")
print(f"  Delta: {new_milw-old_milw:+.2f}pp")

print(f"\nSign convention fix:")
print(f"  OLD column: 'efficacy_lost_to_structural_delay_pp'")
print(f"    Hartford: +0.70 (confusing — positive means 'lost'?)")
print(f"  NEW column: 'efficacy_gain_vs_24h_counterfactual_pp'")
print(f"    Hartford: {df_cost[df_cost['city']=='Hartford']['efficacy_gain_vs_24h_counterfactual_pp'].values[0]:+.2f}pp (negative = worse than 24h baseline)")
print(f"    Milwaukee: {df_cost[df_cost['city']=='Milwaukee']['efficacy_gain_vs_24h_counterfactual_pp'].values[0]:+.2f}pp (positive = better than 24h baseline)")

print(f"\nPerelson anchor check:")
print(f"  MIN_INTEGRATION_COMPLETE = {MIN_INTEGRATION_COMPLETE}h ≈ eclipse phase 0.9d = 21.6h ✓")
print(f"  MIN_SEEDING_MIDPOINT     = {MIN_SEEDING_MIDPOINT}h ≈ virion t½ 0.24d = 5.8h   ✓")
print(f"  REFERENCE_INTEGRATION_COMPLETE = {REFERENCE_INTEGRATION_COMPLETE}h (between eclipse 22h and τ 62h) ✓")
print(f"  VL_UNSUPPRESSED_MEAN_LOG10 = 4.5 → cite NHBS, NOT Perelson ✓")

print("\nAll CSVs written to /mnt/user-data/outputs/")
