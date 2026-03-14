"""
AIDSVu City Profile Parser & City-Specific VL Distribution Generator
=====================================================================

Parses all 34 AIDSVu city datasets and derives empirically-grounded
parameters for StochasticPEPModel, replacing national assumptions
with city-level surveillance data.

THEORETICAL BRIDGE TO PEP MODEL:
    The stochastic PEP model requires two city-specific inputs:

    1. VL distribution parameters (mean_log10, sd_log10):
       Derived from viral suppression rate.
       Logic: If VS% = p, then (1-p) fraction are unsuppressed.
       Unsuppressed VL distribution ~ LogNormal(4.5, 1.1) [Perelson 1996]
       Suppressed VL distribution  ~ LogNormal(1.3, 0.4) [<20-200 copies/mL]
       City mixture mean = weighted average of both components.

    2. Structural delay parameters (hours):
       Derived from late diagnosis % and linkage-to-care %.
       Late dx % → proxy for how long PWID avoid the healthcare system.
       Linkage gap → proxy for infrastructure barrier after diagnosis.
       Combined → city-specific structural delay distribution.

PERELSON 1996 GROUNDING:
    Perelson AS et al. "HIV-1 dynamics in vivo: virion clearance rate,
    infected cell life-span, and viral generation time."
    Science. 1996;271(5255):1582-1586. doi:10.1126/science.271.5255.1582

    Key parameters used here:
    - Viral doubling time in acute infection: ~1.4 days
    - Steady-state VL in untreated chronic infection: 10^3 - 10^6 copies/mL
    - Log-normal VL distribution assumption grounded in this dynamics model
    - Integration timeline (RT → nuclear import → proviral integration):
      6-24h for parenteral, consistent with Perelson's within-host kinetics

Author: AC Demidont, DO, AAHIVS
Date: March 2026
Updated: 2026-03-10 15:32
"""

import pandas as pd
import numpy as np
import glob
import os
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

warnings.filterwarnings('ignore', category=UserWarning)

# --- Path Resolution Utility ---
PROJECT_ROOT = Path(__file__).resolve().parent

def resolve_project_path(relative_path: str) -> str:
    """Resolve a path relative to the project root."""
    return str(PROJECT_ROOT / relative_path)

# =============================================================================
# CONSTANTS: Perelson 1996 + empirical VL distribution parameters
# =============================================================================

# Unsuppressed HIV VL: log-normal distribution
# Source: NHBS/NHAS PWID surveillance; Metzger et al.; Degenhardt et al.
# Median ~30,000 copies/mL → log10 = 4.5 in unsuppressed PWID.
# NOTE: Perelson 1996 (PMID 8599114) grounds the KINETIC PARAMETERS
# (virion t1/2, eclipse phase, generation time) — NOT this VL center.
# Perelson patients had mean VL ~216,000 (log10≈5.3), appropriate for
# kinetics but not PWID community distribution. Cite NHBS for VL=4.5.
VL_UNSUPPRESSED_MEAN_LOG10 = 4.5  # NHBS PWID median, NOT from Perelson
VL_UNSUPPRESSED_SD_LOG10 = 1.1  # Metzger et al. / Degenhardt et al.

# Suppressed HIV VL on effective ART: typically <20-200 copies/mL
# Using 1.3 (=20 copies/mL, LLoD) as effective mean with tight SD
VL_SUPPRESSED_MEAN_LOG10 = 1.5  # ~30 copies/mL effective mean
VL_SUPPRESSED_SD_LOG10 = 0.4  # Very tight: most truly suppressed

# Structural delay mapping: late diagnosis % → hours delay
# Calibration:
#   Late dx <15% → low-barrier city → ~6h structural delay
#   Late dx 15-25% → moderate barrier → ~18h
#   Late dx 25-35% → high barrier → ~30h
#   Late dx >35%  → severe barrier → ~42h
# These map to the pwid_barrier_analysis() scenarios in PEP_parenteral_perelson.py
LATE_DX_DELAY_SLOPE = 1.2  # hours of structural delay per % late diagnosis
LATE_DX_DELAY_INTERCEPT = -12.0  # offset: calibrated to 0h at ~10% late dx

# Linkage gap penalty: each % below 90% linkage adds structural delay
# Rationale: poor linkage = fewer PEP access points = longer navigation time
LINKAGE_GAP_PENALTY_PER_PCT = 0.3  # hours per % below 90%

# IDU transmission columns (City Profile)
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
# PARSER
# =============================================================================

def parse_city_profile(filepath: str) -> Optional[Dict]:
    """
    Extract City Profile variables from a single AIDSVu xlsx file.

    Returns dict of key variables, or None if parse fails.
    """
    city_name = os.path.basename(filepath).split('_AIDSVu')[0]

    try:
        df = pd.read_excel(filepath, sheet_name='City Profile', header=None)

        # Header is at row index 3, data at row index 4
        headers = list(df.iloc[3].values)
        data = list(df.iloc[4].values)

        # Build lookup dict: clean header -> value
        lookup = {}
        for h, v in zip(headers, data):
            if pd.notna(h):
                lookup[str(h)] = v

        def get(col, default=np.nan):
            v = lookup.get(col, default)
            if v == -1 or v == '-1':
                return np.nan  # AIDSVu uses -1 for suppressed/unavailable
            try:
                return float(v)
            except (ValueError, TypeError):
                return default

        def get_str(col, default=''):
            return str(lookup.get(col, default))

        # IDU prevalence: combine male + female weighted by sex distribution
        male_idu_prev = get(COL_MALE_IDU_PREV)
        female_idu_prev = get(COL_FEMALE_IDU_PREV)
        male_pct = get('2023 Male\nPrevalence\nPercent', 75) / 100
        female_pct = 1 - male_pct

        if pd.notna(male_idu_prev) and pd.notna(female_idu_prev):
            idu_prevalence_pct = (male_idu_prev * male_pct +
                                  female_idu_prev * female_pct)
        elif pd.notna(male_idu_prev):
            idu_prevalence_pct = male_idu_prev * male_pct
        else:
            idu_prevalence_pct = np.nan

        # IDU new diagnoses
        male_idu_dx = get(COL_MALE_IDU_DX)
        female_idu_dx = get(COL_FEMALE_IDU_DX)
        if pd.notna(male_idu_dx) and pd.notna(female_idu_dx):
            idu_new_dx_pct = (male_idu_dx * male_pct +
                              female_idu_dx * female_pct)
        else:
            idu_new_dx_pct = np.nan

        return {
            'city': city_name,
            'city_name': get_str('City\nName', city_name),
            'state': get_str(COL_STATE),
            'region': get_str(COL_REGION),
            'population': get(COL_POPULATION),
            'prevalence_cases': get(COL_PREVALENCE_CASES),
            'prevalence_rate': get(COL_PREVALENCE_RATE),
            'new_dx_rate': get(COL_NEW_DX_RATE),
            'new_dx_cases': get(COL_NEW_DX_CASES),
            'viral_suppression_pct': get(COL_VS_PCT),
            'late_dx_pct': get(COL_LATE_DX_PCT),
            'linkage_to_care_pct': get(COL_LINKAGE_PCT),
            'idu_prevalence_pct': idu_prevalence_pct,
            'idu_new_dx_pct': idu_new_dx_pct,
            'male_idu_prev_pct': male_idu_prev,
            'female_idu_prev_pct': female_idu_prev,
            'gini': get(COL_GINI),
            'poverty_pct': get(COL_POVERTY_PCT),
            'housing_burden_pct': get(COL_HOUSING_BURDEN),
        }

    except Exception as e:
        print(f"  WARNING: Failed to parse {city_name}: {e}")
        return None


# =============================================================================
# VL DISTRIBUTION DERIVATION
# =============================================================================

def derive_vl_distribution(vs_pct: float,
                           idu_prevalence_pct: float = np.nan
                           ) -> Dict:
    """
    Derive city-specific VL distribution parameters from viral suppression %.

    Method: Mixture model
        P(VL | city) = VS% * P(VL | suppressed) + (1-VS%) * P(VL | unsuppressed)

    Both components are log-normal (Perelson 1996 within-host dynamics).

    The mixture mean on the log10 scale:
        mu_mix = VS_frac * mu_supp + (1-VS_frac) * mu_unsupp

    The mixture SD (law of total variance):
        var_mix = VS_frac*(sd_supp^2 + mu_supp^2)
                + (1-VS_frac)*(sd_unsupp^2 + mu_unsupp^2)
                - mu_mix^2

    Args:
        vs_pct: Viral suppression percentage (0-100)
        idu_prevalence_pct: IDU-specific prevalence % (modulates unsuppressed
                            component — PWID have lower suppression rates)

    Returns:
        Dict with mean_log10, sd_log10, and interpretation metadata
    """
    vs_frac = vs_pct / 100.0
    unsupp_frac = 1.0 - vs_frac

    # IDU adjustment: PWID typically have ~10-15% lower VS than general PLWH
    # If IDU prevalence data available, apply modest downward correction
    if pd.notna(idu_prevalence_pct) and idu_prevalence_pct > 0:
        # Each % of IDU in prevalence reduces effective suppression by 0.12%
        # Based on NHBS/NHAS data showing PWID VS gaps
        idu_correction = min(idu_prevalence_pct * 0.12, 12.0)  # cap at 12%
        vs_frac_pwid = max(vs_frac - (idu_correction / 100.0), 0.1)
        unsupp_frac_pwid = 1.0 - vs_frac_pwid
    else:
        vs_frac_pwid = vs_frac
        unsupp_frac_pwid = unsupp_frac

    # Mixture mean (log10 scale)
    mean_log10 = (vs_frac_pwid * VL_SUPPRESSED_MEAN_LOG10 +
                  unsupp_frac_pwid * VL_UNSUPPRESSED_MEAN_LOG10)

    # Mixture variance (law of total variance)
    e_x2_supp = VL_SUPPRESSED_SD_LOG10 ** 2 + VL_SUPPRESSED_MEAN_LOG10 ** 2
    e_x2_unsupp = VL_UNSUPPRESSED_SD_LOG10 ** 2 + VL_UNSUPPRESSED_MEAN_LOG10 ** 2
    var_mix = (vs_frac_pwid * e_x2_supp +
               unsupp_frac_pwid * e_x2_unsupp) - mean_log10 ** 2
    sd_log10 = max(np.sqrt(var_mix), 0.3)  # floor at 0.3

    # Classification
    if mean_log10 < 2.5:
        vl_profile = 'low_burden'
    elif mean_log10 < 3.5:
        vl_profile = 'moderate_burden'
    elif mean_log10 < 4.2:
        vl_profile = 'high_burden'
    else:
        vl_profile = 'very_high_burden'

    return {
        'mean_log10': round(mean_log10, 3),
        'sd_log10': round(sd_log10, 3),
        'vs_frac_used': round(vs_frac_pwid, 3),
        'unsupp_frac_used': round(unsupp_frac_pwid, 3),
        'vl_profile': vl_profile,
        'median_vl_copies_ml': round(10 ** mean_log10, 0),
    }


# =============================================================================
# STRUCTURAL DELAY DERIVATION
# =============================================================================

def derive_structural_delay(late_dx_pct: float,
                            linkage_pct: float,
                            gini: float = np.nan,
                            poverty_pct: float = np.nan
                            ) -> Dict:
    """
    Derive city-specific structural delay parameters.

    Components:
    1. Late diagnosis delay: late_dx_pct → hours of avoidance before dx
       Each % above 10% baseline adds LATE_DX_DELAY_SLOPE hours.

    2. Linkage gap penalty: (90 - linkage_pct) → navigation barrier hours
       Each % below 90% adds LINKAGE_GAP_PENALTY_PER_PCT hours.

    3. SDOH modifier: Gini + poverty amplify delay (optional, if available)
       High inequality → more criminalization exposure → longer delay.

    Returns:
        Dict with estimated_delay_hours, delay_category, and components
    """
    # Component 1: Late diagnosis
    late_dx_component = max(
        LATE_DX_DELAY_SLOPE * (late_dx_pct - 10.0) + LATE_DX_DELAY_INTERCEPT,
        0.0
    )

    # Component 2: Linkage gap
    linkage_gap = max(90.0 - linkage_pct, 0.0)
    linkage_component = linkage_gap * LINKAGE_GAP_PENALTY_PER_PCT

    # Component 3: SDOH amplifier (0-6h bonus)
    sdoh_component = 0.0
    if pd.notna(gini) and pd.notna(poverty_pct):
        # Gini > 0.45 = high inequality; poverty > 20% = high poverty
        gini_factor = max(gini - 0.40, 0) * 20.0  # 0-4h
        poverty_factor = max(poverty_pct - 15.0, 0) * 0.1  # 0-2h
        sdoh_component = min(gini_factor + poverty_factor, 6.0)

    total_delay = late_dx_component + linkage_component + sdoh_component

    # Category
    if total_delay < 8:
        delay_category = 'low_barrier'
    elif total_delay < 18:
        delay_category = 'moderate_barrier'
    elif total_delay < 30:
        delay_category = 'high_barrier'
    else:
        delay_category = 'severe_barrier'

    return {
        'estimated_structural_delay_h': round(total_delay, 1),
        'delay_category': delay_category,
        'late_dx_component_h': round(late_dx_component, 1),
        'linkage_component_h': round(linkage_component, 1),
        'sdoh_component_h': round(sdoh_component, 1),
    }


# =============================================================================
# MASTER PARSER: All 34 cities
# =============================================================================

def build_city_profiles(data_dir: Optional[str] = None) -> pd.DataFrame:
    """
    Parse all AIDSVu files and build master city profile DataFrame.
    """
    if data_dir is None:
        data_dir = resolve_project_path('aidsvu datasets')

    files = sorted(glob.glob(os.path.join(data_dir, '*AIDSVu*.xlsx')))
    print(f"Found {len(files)} AIDSVu files. Parsing...")
    rows = []
    for fp in files:
        city_name = os.path.basename(fp).split('_AIDSVu')[0]
        print(f"  Parsing {city_name}...", end=' ')

        profile = parse_city_profile(fp)
        if profile is None:
            print("FAILED")
            continue

        # Derive VL distribution
        vl_params = derive_vl_distribution(
            vs_pct=profile['viral_suppression_pct'],
            idu_prevalence_pct=profile['idu_prevalence_pct']
        )

        # Derive structural delay
        delay_params = derive_structural_delay(
            late_dx_pct=profile['late_dx_pct'],
            linkage_pct=profile['linkage_to_care_pct'],
            gini=profile['gini'],
            poverty_pct=profile['poverty_pct']
        )

        # Expected annual PWID PEP-eligible exposures
        # Rough estimate: new_dx_cases * idu_new_dx_fraction * exposure_multiplier
        # Exposure multiplier ~5x: most PWID exposures don't result in dx
        idu_dx_frac = (profile['idu_new_dx_pct'] or 0) / 100.0
        estimated_pwid_exposures = (
                (profile['new_dx_cases'] or 0) * idu_dx_frac * 5.0
        )

        row = {**profile, **vl_params, **delay_params,
               'estimated_annual_pwid_exposures': round(estimated_pwid_exposures) if pd.notna(
                   estimated_pwid_exposures) else 0}
        rows.append(row)
        print(f"VS={profile['viral_suppression_pct']:.0f}% | "
              f"LateDx={profile['late_dx_pct']:.0f}% | "
              f"IDU={profile['idu_prevalence_pct']:.1f}% | "
              f"VL_mean_log10={vl_params['mean_log10']:.2f} | "
              f"Delay={delay_params['estimated_structural_delay_h']:.0f}h")

    df = pd.DataFrame(rows)

    return df


# =============================================================================
# CITY-STRATIFIED STOCHASTIC PEP ANALYSIS
# =============================================================================

def run_city_stratified_analysis(city_df: pd.DataFrame,
                                 pep_delay_hours: float = 24.0,
                                 n_simulations: int = 5000
                                 ) -> pd.DataFrame:
    """
    Run stochastic PEP model for each city using empirically derived parameters.

    For each city:
        1. Uses city-specific VL distribution (from VS%)
        2. Uses city-specific structural delay (from late dx% + linkage%)
        3. Computes expected PEP efficacy at that city's structural delay
        4. Computes P(futile) — proportion where PEP is biologically too late

    Args:
        city_df: Output of build_city_profiles()
        pep_delay_hours: Override delay (use city-specific if None)
        n_simulations: Monte Carlo samples per city

    Returns:
        DataFrame with city-level PEP efficacy estimates
    """
    # Import here to avoid circular dependency
    from PEP_parenteral_perelson import ParenteralExposureModel

    rng = np.random.default_rng(42)
    results = []

    print(f"\nRunning city-stratified analysis "
          f"(n={n_simulations} per city)...")

    for _, city in city_df.iterrows():
        city_name = city['city']
        mean_log10 = city['mean_log10']
        sd_log10 = city['sd_log10']

        # Use city structural delay unless override provided
        delay_h = (city['estimated_structural_delay_h']
                   if pep_delay_hours is None
                   else pep_delay_hours)

        # Sample VL from city distribution
        sampled_log10 = rng.normal(loc=mean_log10, scale=sd_log10,
                                   size=n_simulations)
        sampled_log10 = np.clip(sampled_log10, 1.3, 7.0)
        sampled_vl = 10 ** sampled_log10

        # Compute efficacy for each sampled VL
        efficacies = []
        for vl in sampled_vl:
            model = ParenteralExposureModel(source_viral_load=vl)
            eff = model.pep_efficacy(delay_h)['pep_efficacy']
            efficacies.append(eff)

        efficacies = np.array(efficacies)

        results.append({
            'city': city_name,
            'state': city['state'],
            'region': city['region'],
            'viral_suppression_pct': city['viral_suppression_pct'],
            'late_dx_pct': city['late_dx_pct'],
            'linkage_to_care_pct': city['linkage_to_care_pct'],
            'idu_prevalence_pct': city['idu_prevalence_pct'],
            'vl_mean_log10': mean_log10,
            'vl_sd_log10': sd_log10,
            'vl_profile': city['vl_profile'],
            'structural_delay_h': delay_h,
            'delay_category': city['delay_category'],
            'pep_mean_efficacy_pct': round(np.mean(efficacies) * 100, 2),
            'pep_median_efficacy_pct': round(np.median(efficacies) * 100, 2),
            'pep_ci95_lower_pct': round(np.percentile(efficacies, 2.5) * 100, 2),
            'pep_ci95_upper_pct': round(np.percentile(efficacies, 97.5) * 100, 2),
            'p_futile_pct': round(np.mean(efficacies < 0.1) * 100, 2),
            'p_below_50pct': round(np.mean(efficacies < 0.5) * 100, 2),
            'estimated_pwid_exposures': city['estimated_annual_pwid_exposures'],
            'new_dx_cases': city['new_dx_cases'],
            'prevalence_rate': city['prevalence_rate'],
            'gini': city['gini'],
            'poverty_pct': city['poverty_pct'],
        })

        print(f"  {city_name:<18} | delay={delay_h:.0f}h | "
              f"E[eff]={np.mean(efficacies) * 100:.1f}% | "
              f"P(futile)={np.mean(efficacies < 0.1) * 100:.1f}%")

    return pd.DataFrame(results).sort_values('pep_mean_efficacy_pct',
                                             ascending=True)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("AIDSVu CITY PROFILE PARSER")
    print("Prevention Theorem — City-Stratified Stochastic PEP Model")
    print("=" * 70)

    output_dir = resolve_project_path('results/city_analysis')
    os.makedirs(output_dir, exist_ok=True)

    # --- Parse all cities ---
    city_df = build_city_profiles()

    # Define cleaning function for export only
    def clean_df_for_export(df):
        # Create a copy to avoid modifying the original dataframe
        df_clean = df.copy()
        # Clean the cell values
        df_clean = df_clean.replace(to_replace=r'\n|\r', value=' ', regex=True)
        # Clean the index if it's an object/string type (CRITICAL for San Juan / New Haven)
        if df_clean.index.dtype == 'object':
            df_clean.index = df_clean.index.str.replace(r'\n|\r', ' ', regex=True)
        return df_clean

    # Save master profile
    path1 = os.path.join(output_dir, 'city_vl_profiles.csv')
    clean_df_for_export(city_df).to_csv(path1, index=False)
    print(f"\nSaved master city profiles: {path1} ({len(city_df)} cities)")

    # --- Summary statistics ---
    print("\n" + "=" * 70)
    print("CITY PROFILE SUMMARY")
    print("=" * 70)
    # Use clean_df_for_export for summary to avoid hidden newlines in printed tables
    # but keep city_df stable for calculations
    city_df_summary = clean_df_for_export(city_df)

    print(f"\nViral suppression range: "
          f"{city_df['viral_suppression_pct'].min():.1f}% – "
          f"{city_df['viral_suppression_pct'].max():.1f}%")
    print(f"Late diagnosis range:     "
          f"{city_df['late_dx_pct'].min():.1f}% – "
          f"{city_df['late_dx_pct'].max():.1f}%")
    print(f"IDU prevalence range:     "
          f"{city_df['idu_prevalence_pct'].min():.1f}% – "
          f"{city_df['idu_prevalence_pct'].max():.1f}%")
    print(f"VL mean_log10 range:      "
          f"{city_df['mean_log10'].min():.2f} – "
          f"{city_df['mean_log10'].max():.2f}")
    print(f"Structural delay range:   "
          f"{city_df['estimated_structural_delay_h'].min():.0f}h – "
          f"{city_df['estimated_structural_delay_h'].max():.0f}h")

    print("\nDelay category distribution:")
    print(city_df_summary['delay_category'].value_counts().to_string())
    print("\nVL profile distribution:")
    print(city_df_summary['vl_profile'].value_counts().to_string())

    # Top 5 highest burden cities
    print("\nTop 5 cities by structural delay (worst barriers):")
    worst = city_df_summary.nlargest(5, 'estimated_structural_delay_h')[
        ['city', 'state', 'late_dx_pct', 'linkage_to_care_pct',
         'estimated_structural_delay_h', 'delay_category', 'mean_log10']]
    print(worst.to_string(index=False))

    print("\nTop 5 cities by VL burden (highest community VL):")
    highest_vl = city_df_summary.nlargest(5, 'mean_log10')[
        ['city', 'state', 'viral_suppression_pct', 'idu_prevalence_pct',
         'mean_log10', 'vl_profile']]
    print(highest_vl.to_string(index=False))

    # --- City-stratified PEP analysis ---
    # Run at 24h (optimal) and city-specific structural delay
    print("\n" + "=" * 70)
    print("CITY-STRATIFIED STOCHASTIC PEP ANALYSIS")
    print("Using city-specific structural delay as PEP initiation time")
    print("=" * 70)

    pep_results = run_city_stratified_analysis(
        city_df=city_df,
        pep_delay_hours=None,  # Use city-specific structural delay
        n_simulations=5000
    )

    path2 = os.path.join(output_dir, 'city_pep_efficacy_results.csv')
    clean_df_for_export(pep_results).to_csv(path2, index=False)
    print(f"\nSaved city PEP results: {path2}")

    # Print ranked results
    print("\n" + "=" * 70)
    print("CITY-LEVEL EXPECTED PEP EFFICACY (ranked worst → best)")
    print("=" * 70)
    cols = ['city', 'state', 'structural_delay_h', 'pep_mean_efficacy_pct',
            'p_futile_pct', 'vl_profile', 'delay_category']
    print(pep_results[cols].to_string(index=False))

    # --- Compare at 24h to show structural delay impact ---
    print("\n" + "=" * 70)
    print("COUNTERFACTUAL: All cities at 24h delay (no structural barrier)")
    print("=" * 70)
    pep_24h = run_city_stratified_analysis(
        city_df=city_df,
        pep_delay_hours=24.0,
        n_simulations=5000
    )
    path3 = os.path.join(output_dir, 'city_pep_efficacy_24h_counterfactual.csv')
    clean_df_for_export(pep_24h).to_csv(path3, index=False)

    # Merge to show structural delay cost
    merged = pep_results[['city', 'pep_mean_efficacy_pct', 'structural_delay_h']].merge(
        pep_24h[['city', 'pep_mean_efficacy_pct']].rename(
            columns={'pep_mean_efficacy_pct': 'pep_24h_counterfactual_pct'}),
        on='city'
    )
    merged['efficacy_lost_to_structural_delay_pp'] = (
            merged['pep_24h_counterfactual_pct'] - merged['pep_mean_efficacy_pct']
    ).round(2)
    merged = merged.sort_values('efficacy_lost_to_structural_delay_pp',
                                ascending=False)

    path4 = os.path.join(output_dir, 'city_structural_delay_cost.csv')
    clean_df_for_export(merged).to_csv(path4, index=False)

    print("\nEfficacy lost due to structural delay (worst cities):")
    print(merged.head(10).to_string(index=False))

    print(f"\n{'=' * 70}")
    print("OUTPUT FILES:")
    print(f"  {path1}")
    print(f"  {path2}")
    print(f"  {path3}")
    print(f"  {path4}")
    print(f"{'=' * 70}")
    print("""
PERELSON 1996 GROUNDING NOTE:
    VL distribution parameters are grounded in Perelson et al. (1996)
    within-host HIV dynamics model. The log-normal VL distribution with
    mean log10=4.5 for untreated infection reflects the paper's estimate
    of steady-state viral burden during chronic infection.

    The integration timeline compression (parenteral: 4-60h vs mucosal:
    72-120h) is consistent with Perelson's viral kinetics: at high inoculum,
    reverse transcription initiates within hours of entry, compressing the
    window within which PEP can achieve R(0)=0.

    Citation: Perelson AS, Neumann AU, Markowitz M, Leonard JM, Ho DD.
    HIV-1 dynamics in vivo: virion clearance rate, infected cell life-span,
    and viral generation time. Science. 1996;271(5255):1582-6.
    doi: 10.1126/science.271.5255.1582
    """)


if __name__ == "__main__":
    main()