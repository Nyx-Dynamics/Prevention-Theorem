# Finite Prevention Windows for HIV Post-Exposure Prophylaxis
## Science Advances Submission — Manuscript ID: aeh5879 (v2)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20044747.svg)](https://doi.org/10.5281/zenodo.20044747)

**Title:** Finite Prevention Windows for HIV Post-Exposure Prophylaxis: Irreversible Proviral Integration Defines Route-Specific Intervention Limits  
**Author:** A.C. Demidont, Nyx Dynamics LLC  
**Date prepared:** 2026-05-05  
**Zenodo:** https://doi.org/10.5281/zenodo.20044747 (v2-submission placeholder)

---

## Quickstart: Reproducing Results

To reproduce all figures and numerical claims for the v2 submission:

1. **Environment Setup** (Python 3.9+ recommended):
   ```bash
   pip install -r requirements.txt
   ```

2. **Execute Reproduction Pipeline**:
   ```bash
   python3 reproduce_all_v2.py
   ```
   *Note: This script sequentially generates Figures 1-4 and their underlying data CSVs, saving them to a unique, timestamped directory under `v2_revision/runs/`.*

3. **Using Makefile** (Alternative):
   ```bash
   make all
   ```

---

## File Inventory (v2)

| File / Directory | Purpose |
|---|---|
| `reproduce_all_v2.py` | Master script for end-to-end reproduction. |
| `v2_revision/` | Primary directory for v2 submission artifacts. |
| `v2_revision/numerical_claims_v2.csv` | Reconciled summary of all numerical claims. |
| `SRC/multiscale_model/` | Within-host multiscale model (Phase 1 & 2). |
| `SRC/perelson/stochastic/` | Stochastic PEP efficacy & VL uncertainty analysis. |
| `v2_revision/city_stratified_figures.py` | City-stratified analysis for 34 US cities. |
| `archive/` | Stale scripts, data, and documentation from v1. |

---

## Order of Operations (v2 Workflow)

1. **Model Execution**: Figures are generated from either cached results in `SRC/multiscale_model/results_v3/` or via fresh Monte Carlo simulations.
2. **Result Organization**: `reproduce_all_v2.py` captures all PNG/PDF/CSV outputs and organizes them into timestamped run folders.
3. **Verification**: Numerical outputs (e.g., $t_{\text{crit}}$ of 34.5h parenteral / 60.5h mucosal) are printed to the console and compared against `v2_revision/numerical_claims_v2.csv`.

---

## Software Environment

This repository is pinned for high-fidelity reproduction. 
- **Python Version**: 3.9.13
- **Primary Dependencies**: `numpy`, `scipy`, `matplotlib`, `pandas`, `seaborn`, `openpyxl`.
- **Reproducibility**: Uses deterministic seeds and a non-interactive Matplotlib backend (`Agg`).

## Tags and Editorial Transparency

- `disclosure-2026-04-30` (commit `d047d2d`) — preserved as the editorial transparency snapshot at the time of the self-disclosed correction (October 30, 2026 disclosure).
- `v2-prep-2026-05-05` (commit `37e27ea`) — the SHA cited in the v2 manuscript's Methods and Data-Availability sections; numerical claims reproduce from this commit.
- `v2-submission-2026-05-05` — the polished submission state with filled DOIs, comprehensive source-code inclusion, and the LaTeX deliverables in `v2_revision/`.

---

## License
This project is licensed under the terms of the license included in the repository.
