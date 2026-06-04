# The HIV Post-Exposure Prophylaxis Window: A Multiscale Framework Linking Within-Host Integration Kinetics to Population-Scale Structural Access
## PLOS Computational Biology submission

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20044747.svg)](https://doi.org/10.5281/zenodo.20044747)

**Author:** A.C. Demidont, DO, AAHIVS — Nyx Dynamics LLC  
**Zenodo (v2 baseline):** https://doi.org/10.5281/zenodo.20044747  
**Zenodo (v3):** *DOI pending — will be minted upon submission confirmation*

---

## Quickstart: Reproducing Results

**v3 full pipeline** (multiscale baseline + PK-driven framework + pharmacy sensitivity):

```bash
pip install -r requirements.txt
python3 reproduce_all_v3.py
```

**v2 baseline only** (canonical multiscale baseline, unchanged from v2):

```bash
python3 reproduce_all_v2.py
# or: make all
```

Both scripts output to a timestamped run directory. The v3 script produces all main-text figures plus S14 pharmacy sensitivity outputs. v2 outputs remain the canonical reference for Figures 1–3 and S1.

---

## File Inventory

| File / Directory | Purpose |
|---|---|
| `reproduce_all_v3.py` | Master script: v2 baseline + v3 PK-driven layer + pharmacy sensitivity |
| `reproduce_all_v2.py` | v2-only master script; canonical multiscale baseline reproduction |
| `v3_revision/` | v3 submission artifacts |
| `v3_revision/r3_pk_pd/` | PK-driven framework: `pk_model.py`, `pd_model.py`, `effective_epsilon.py`, `pharmacy_sensitivity.py`, `test_regression.py` |
| `v3_revision/data/Table_34_cities_full.csv` | 34-city AIDSVu-derived input for pharmacy sensitivity (S14) |
| `v3_revision/numerical_claims_v3.csv` | All verifiable numerical claims for v3 (v2 baseline + PK/pharmacy additions) |
| `v2_revision/` | v2 submission artifacts (Science Advances aeh5879 — closed); preserved for disclosure trail |
| `v2_revision/numerical_claims_v2.csv` | v2 reconciled numerical claims |
| `SRC/multiscale_model/` | Within-host multiscale model (Phase 1 & 2); provides t_crit and VL realizations |
| `SRC/perelson/stochastic/` | Stochastic PEP efficacy and VL uncertainty analysis |
| `stochastic_layers/` | Stochastic envelope analyses |
| `route_models/` | Route-dependent compression scaffolding |
| `core_theorem/` | Mathematical core (theorems and proofs) |
| `aidsvu datasets/` | AIDSVu county/state surveillance data (34 cities + state viral suppression) |
| `results/` | Output directory for city analysis and integrated results |
| `archive_exploratory/` | Stale v1/exploratory artifacts |
| `plos_comp_bio/` | PLOS Comp Bio submission package: manuscript TeX, supplement, figures |
| `audit_check.py` | v2 reproducibility audit |
| `audit_check_v3.py` | v3 reproducibility audit (PK-driven + pharmacy sensitivity claims) |
| `Makefile` | `make all` reproduces v2 baseline |

---

## v3 Contribution Layer: Order of Operations

The v3 pipeline is two-part. v3 **does not replace** v2 — it adds the PK-driven framework as a calibration and sensitivity layer on top of the unchanged multiscale baseline.

**Part 1 — v2 baseline reproduction:**  
Calls existing v2 scripts to regenerate Figures 1, 2, 3, S1. v3 main-text figures are identical to v2 at perfect adherence because the PK-driven framework recovers the multiscale baseline by design.

**Part 2 — v3-specific analyses:**
1. `v3_revision/r3_pk_pd/test_regression.py` — validates that the PK-driven framework recovers the multiscale t_crit and envelope within calibration tolerance at ρ=1.0
2. `v3_revision/r3_pk_pd/pharmacy_sensitivity.py` — generates Supplement Figure S14: pharmacy access delay sweep (Δt ∈ {0,2,4,6,8,10,12} h) across 34-city panel

All outputs land in `v3_revision/runs/run_<timestamp>/`.

---

## Numerical Claims and Verification

Core values this pipeline must reproduce (tolerances in parentheses):

| Claim | Value | Tolerance |
|---|---|---|
| Parenteral t_crit (η=0.05, V₀=10³, CV=0.3) | 34.5 h | ±0.5 h |
| Mucosal t_crit (η=0.05, V₀=1, CV=0.3) | 60.5 h | ±0.5 h |
| Compression ratio (mucosal/parenteral) | ~1.75× | — |
| Envelope bound (median=96h, GSD=2.0) | ~11.3% | ±0.5 pp |
| PK-driven t_crit at ρ=1.0 (parenteral) | 34.0 h | ±0.5 h |
| PK-driven t_crit at ρ=1.0 (mucosal) | 59.5 h | ±0.5 h |
| PK-driven envelope at ρ=1.0 | 11.2% | ±0.5 pp |
| Hartford E_PEP at Δt_pharmacy = 8 h | < 0.05 | — |
| Envelope across pharmacy sweep | 11.5% → 9.6% | monotone ↓ |
| Mean E_PEP across 34 cities (Δt=0 → 12h) | 96.1% → 83.7% | — |

Full claim table: `v3_revision/numerical_claims_v3.csv`.  
Run `python3 audit_check_v3.py` to verify all claims against current outputs.

---

## Software Environment

- **Python Version**: 3.9.13
- **Primary Dependencies**: `numpy`, `scipy`, `matplotlib`, `pandas`, `seaborn`, `openpyxl`
- **Reproducibility**: Deterministic seeds; non-interactive Matplotlib backend (`Agg`)

---

## Tags and Editorial Transparency

- `disclosure-2026-04-30` (commit `d047d2d`) — editorial transparency snapshot at time of voluntary numerical-error disclosure
- `v1.3.0` — Science submission v5 (aeh1546)
- `v2-prep-2026-05-05` (commit `37e27ea`) — SHA cited in v2 manuscript Methods; multiscale baseline numerical claims reproduce from this commit
- `v2-submission-2026-05-05` — polished Science Advances submission state (closed submission)
- `v3-prep-2026-06-03` — SHA to be cited in v3 manuscript Methods upon PR merge *(pending AC confirmation)*
- `v3-submission-YYYY-MM-DD` — will be tagged at PLOS Comp Bio submission *(pending)*

---

## License

MIT — see `LICENSE`.
