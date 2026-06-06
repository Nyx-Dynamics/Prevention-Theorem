# The HIV Post-Exposure Prophylaxis Window: A Multiscale Framework Linking Within-Host Integration Kinetics to Population-Scale Structural Access

**PLOS Computational Biology** — submission PCOMPBIOL-S-26-01758-2 (v4)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20559022.svg)](https://doi.org/10.5281/zenodo.20559022)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20044747.svg)](https://doi.org/10.5281/zenodo.20044747)

**Author:** A.C. Demidont, DO, AAHIVS — Nyx Dynamics LLC  
**Zenodo (v4, current):** https://doi.org/10.5281/zenodo.20559022  
**Zenodo (v2 baseline):** https://doi.org/10.5281/zenodo.20044747

---

## Quickstart: Reproducing Results

**v4 full pipeline** (multiscale baseline + PK-driven framework + envelope-corridor pharmacy analysis):

```bash
pip install -r requirements.txt
python3 reproduce_all_v4.py
```

**v2 baseline only** (canonical multiscale baseline, unchanged from v2):

```bash
python3 reproduce_all_v2.py
# or: make all
```

Both scripts write all outputs to a timestamped run directory under `v3_revision/runs/`. v4 also writes a `run_metadata.txt` provenance file to each run.

---

## File Inventory

| File / Directory | Purpose |
|---|---|
| `reproduce_all_v4.py` | Master v4 pipeline: multiscale baseline + PK-driven framework + envelope-corridor pharmacy |
| `reproduce_all_v2.py` | v2-only master; canonical multiscale baseline |
| `finite_windows_v4/` | Self-contained v4 submission package (see below) |
| `finite_windows_main_plos_comp_bio.tex` | Root-level manuscript (kept in sync with v4 final) |
| `finite_windows_supplement_plos_comp_bio.tex` | Root-level supplement (kept in sync with v4 final) |
| `v3_revision/` | v3/v4 submission artifacts |
| `v3_revision/r3_pk_pd/` | PK-driven framework + v4 corridor scripts (see below) |
| `v3_revision/results/pharmacy_sensitivity_corrected/` | Figure S14 + pharmacy corridor CSVs |
| `v3_revision/runs/` | Timestamped reproduction run archives |
| `v3_revision/numerical_claims_v3.csv` | v3 claims registry |
| `v3_revision/numerical_claims_v4.csv` | v4 claims registry (corridor replaces scalar envelope bound) |
| `v2_revision/` | v2 artifacts (Science Advances aeh5879 — closed); preserved for disclosure trail |
| `SRC/multiscale_model/` | Within-host multiscale model (Phase 1 & 2); provides t_crit and VL realizations |
| `SRC/perelson/stochastic/` | Stochastic PEP efficacy and VL uncertainty analysis |
| `aidsvu datasets/` | AIDSVu county/state surveillance data (34 cities + state viral suppression) |
| `results/` | Output directory for city analysis and integrated results |
| `audit_check_v3.py` | Verifies all v3/v4 numerical claims against current outputs |
| `Makefile` | `make all` reproduces v2 baseline |

### `v3_revision/r3_pk_pd/` scripts

| Script | Purpose |
|---|---|
| `pk_model.py` | PK model: plasma concentration under adherence parameter ρ |
| `pd_model.py` | PD model: stage-dependent drug efficacy ε(t) |
| `effective_epsilon.py` | Combines PK+PD to produce E_PEP corridor curves |
| `pharmacy_sensitivity.py` | Computes envelope corridor + city positions + displacement counts |
| `envelope_corridor_figure.py` | Figure S14: 2D corridor with 34-city overlay and Hartford displacement arrow |
| `generate_S2_tables.py` | Supplementary S2 table generation |
| `test_regression.py` | Verifies PK-driven framework recovers multiscale baseline at ρ=1.0 |

### `finite_windows_v4/` submission package

Self-contained directory mirroring the PLOS Comp Bio submission:

- `finite_windows_v4_main_plos_comp_bio.{tex,pdf}` — manuscript
- `finite_windows_v4_supplement_plos_comp_bio.{tex,pdf}` — supplement
- `PCOMPBIOL-S-26-01758-2_submission.pdf` — submitted PDF bundle
- `Demidont_PLOS_CompBio_cover_letter.docx` — cover letter
- All 5 figures (Figures 1–3, S1, S14)
- All source data CSVs (city analysis, stochastic VL, pharmacy corridor, multiscale)
- All r3_pk_pd scripts + `reproduce_all_v4.py`
- `numerical_claims_v3.csv`, `numerical_claims_v4.csv`
- Archive zips: `finite_windows_zenodo_v4.zip`, `finite_windows_v4_{code,figures}.zip`, `finite_windows_plos_comp_bio_data_numerical_claims.zip`

---

## v4 Correction: Envelope-Corridor Framing

v4 corrects a methodological error in the v3 pharmacy sensitivity output. The prior "envelope bound" scalar conflated three quantities that must remain separate:

- **F_access(t_crit)** — population-level access distribution
- **ε_max** — treated as constant 0.98 rather than a function of remaining window
- **ε_min = 0.05** — a phantom residual floor with no causal basis

Under v4, pharmacy access is modeled as an **upstream gating event**: patients who cannot acquire medication before t_crit receive no drug (E_PEP = 0, not ε_min). The "envelope" is a **2D corridor** in (t_acq, E_PEP) space — bounded above by perfect-adherence PK (ρ=1.0) and below by low-adherence PK (ρ=0.30). Pharmacy delay displaces cities *rightward along a fixed kinetic boundary*; it does not compress the corridor.

**Superseded v3 outputs (do not cite):**
- Scalar `envelope_sweep_dt0 = 11.5%` — Bernoulli mixture; misleads
- Scalar `envelope_sweep_dt12h = 9.6%` — same
- Phantom ε_min = 0.05 floor — no drug → no effect → E_PEP = 0

---

## Pipeline: Order of Operations

v4 is a three-part pipeline. v2 baseline outputs are **unchanged**.

**Part 1 — v2 baseline (Figures 1–3, S1):**  
Calls v2 scripts to regenerate main-text figures. The PK-driven framework recovers the multiscale baseline at ρ=1.0 by design.

**Part 2 — PK-driven framework + v4 corridor pharmacy:**
1. `test_regression.py` — validates PK-driven framework recovers multiscale t_crit within tolerance at ρ=1.0
2. `pharmacy_sensitivity.py` — computes corridor positions + displacement counts for all 34 cities across Δt_pharm ∈ {0,2,4,6,8,10,12} h
3. `envelope_corridor_figure.py` — generates Figure S14

**Part 3 — Numerical claims verification:**  
Prints all v4 claims; checks against `v3_revision/numerical_claims_v4.csv`.

All outputs land in `v3_revision/runs/run_<timestamp>_v4/`.

---

## Numerical Claims (v4)

| Claim | Value | Tolerance |
|---|---|---|
| Parenteral t_crit (η=0.05, V₀=10³, CV=0.3) | 34.5 h | ±0.5 h |
| Mucosal t_crit (η=0.05, V₀=1, CV=0.3) | 60.5 h | ±0.5 h |
| Compression ratio (mucosal/parenteral) | ~1.75× | — |
| F_access × t_crit bound (canonical params) | ~11.3% | ±0.5 pp |
| PK t_crit at ρ=1.0 (parenteral, upper corridor) | 34.0 h | ±0.5 h |
| PK t_crit at ρ=0.30 (parenteral, lower corridor) | 32.0 h | ±0.5 h |
| Hartford at Δt_pharm=0: t_acq, E_PEP_upper | 27.2 h, 0.47 | on corridor |
| Hartford at Δt_pharm=8: t_acq vs t_crit | 35.2 h > 34.0 h | off corridor, E_PEP=0 |
| Cities past t_crit at Δt_pharm ∈ {0,2,4,6} h | 0 of 34 | — |
| Cities past t_crit at Δt_pharm ∈ {8,10,12} h | 1 of 34 (Hartford) | — |
| Cohort mean E_PEP_upper at Δt_pharm=0 | 0.961 | high-burden metros only |
| Cohort mean E_PEP_upper at Δt_pharm=12 h | 0.837 | high-burden metros only |

Full registry: `v3_revision/numerical_claims_v4.csv`.

---

## Software Environment

- **Python:** 3.9.13
- **Dependencies:** `numpy`, `scipy`, `matplotlib`, `pandas`, `seaborn`, `openpyxl`
- **Reproducibility:** Deterministic seeds; non-interactive Matplotlib backend (`Agg`)

---

## Version History and Tags

| Tag / Branch | Commit | Description |
|---|---|---|
| `disclosure-2026-04-30` | `d047d2d` | Editorial transparency snapshot at time of voluntary numerical-error disclosure |
| `v2-prep-2026-05-05` | `37e27ea` | SHA cited in v2 manuscript Methods; multiscale baseline reproduces from here |
| `v2-submission-2026-05-05` | — | Science Advances submission state (closed) |
| `v3-prep-2026-06-03` | current | v4 pipeline + submission package + figure layout fixes |

**Zenodo deposits:**
- v4 (current): https://doi.org/10.5281/zenodo.20559022
- v2 baseline: https://doi.org/10.5281/zenodo.20044747

---

## License

MIT — see `LICENSE`.
