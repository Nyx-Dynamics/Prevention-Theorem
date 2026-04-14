# Finite Prevention Windows for HIV Post-Exposure Prophylaxis
## Science Advances Submission — Manuscript ID: aeh1546

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19011771.svg)](https://doi.org/10.5281/zenodo.19011771)

**Submitted:** March 13, 2026  
**Journal:** *Science Advances*  
**Author:** A.C. Demidont, DO — Nyx Dynamics, LLC  
**Zenodo:** https://doi.org/10.5281/zenodo.19011771

## Contents

### Submission PDF
- `finite_windows_final_Science_combined.pdf` — Complete submission
  package: manuscript, supplement, figures, tables, and figure legends

### Makefile
- `Makefile` — Reproduces all figures and datasets. Run `make all`
  to regenerate. Outputs are timestamped under `outputs/`.

### Reproducible Outputs — `2026-03-13_1941/`
Final `make all` run as submitted to *Science*. Timestamped outputs include:

**Figures**
- `FIgure_1_pep_parenteral_vl_sweep.png` — Fig 1: Parenteral VL sweep
- `Figure_2_pep_stochastic_vl_uncertainty.png` — Fig 2: Stochastic VL uncertainty
- `results/city_analysis/Fig_CityStratified_PEP.png` — Fig 3: City-stratified PEP efficacy
- `results/city_analysis/Fig_CityComparison_Focus.png` — Supp Fig S3: City comparison

**City-level Data (AIDSVu)**
- `city_vl_profiles.csv` — Viral load profiles for 34 US cities
- `city_pep_efficacy_results.csv` — PEP efficacy estimates by city
- `city_pep_efficacy_24h_counterfactual.csv` — 24h access counterfactual
- `city_structural_delay_cost.csv` — Efficacy lost to structural delay

**Stochastic Analysis (`pep_stochastic_results/`)**
- `pep_stochastic_efficacy_curves.csv` — Full efficacy curves
- `pep_stochastic_summary_timepoints.csv` — Summary at key timepoints
- `pep_stochastic_ci_width_regimes.csv` — CI width regime analysis
- `pep_vl_knowledge_premium.csv` — VL knowledge premium analysis

## Reproducing Results
```bash
pip install -r requirements.txt
make all
```
Outputs will be saved to `outputs/<timestamp>/`.
