# aeh5879 v2 — Submission Package

**Manuscript ID:** aeh5879 (Science Advances, transferred from Science)
**Title:** Finite Prevention Windows for HIV Post-Exposure Prophylaxis: Irreversible Proviral Integration Defines Route-Specific Intervention Limits
**Author:** A.C. Demidont, Nyx Dynamics LLC
**Date prepared:** 2026-05-05

---

## File inventory

| File | Purpose | Status |
|---|---|---|
| `main_v2.tex` | Manuscript v2 (LaTeX, Science Advances format) | ready for review |
| `supplementary_v2.tex` | Supplementary Materials v2 (LaTeX) | ready for review |
| `cover_letter_v2.tex` | Cover letter to Erin LaFlame (LaTeX) | ready for review |
| `v2_text_deltas.md` | Section-by-section diff of changes from v1 | reference |
| `numerical_claims_v2.csv` | Reconciled summary of all numerical claims with source files and SHAs | ready (24 rows) |
| `manuscript_v1_text.txt` | Extracted text of v1 manuscript for reference | reference |
| `supplement_v1_text.txt` | Extracted text of v1 supplement for reference | reference |
| `README_v2.md` | This file | ready |

## Compile instructions

```
pdflatex main_v2 && pdflatex main_v2 && pdflatex main_v2
pdflatex supplementary_v2 && pdflatex supplementary_v2
pdflatex cover_letter_v2
```

Three passes ensure cross-references and bibliography resolve. Figures are referenced from a `figures/` subdirectory; copy the v2-regenerated figure PNGs there before compiling.

## Placeholders to fill before submission

The following placeholders appear in the LaTeX sources and **must** be replaced with concrete values before the package goes to Erin:

| Placeholder | What to replace with | When |
|---|---|---|
| `[V2_SUBMISSION_SHA]` | The new annotated tag SHA (e.g., output of `git rev-parse v2-submission-2026-05-XX`) | After v2 commits land and tag is created |
| `[V2_ZENODO_CODE_DOI]` | Code Zenodo DOI from a fresh deposit at the v2-submission tag (NOT the disclosure-era DOI) | After fresh Zenodo deposit |
| `[V2_ZENODO_DATA_DOI]` | Dataset Zenodo DOI from the same fresh deposit | After fresh Zenodo deposit |

The placeholders are intentional — they document at what step in the workflow each value gets pinned. Do not pre-fill any of these from prior values; the v2 deposit and v2 SHA are new artifacts.

## Order of operations (pre-submission checklist)

1. **Read** `main_v2.tex`, `supplementary_v2.tex`, `cover_letter_v2.tex`, and `numerical_claims_v2.csv` end to end. Mark up any wording or numbers to change.
2. **Apply repository hygiene**: pin `requirements.txt` versions, restore `SRC/multiscale_model/run_mc_v3.py` to working tree, verify `README.md` documents file layout and reproduction sequence.
3. **N-sensitivity sweep**: run the model at $V_0=10^3$, CV=0.3, $N \in \{200, 500, 1000, 5000\}$ to confirm $t_{\text{crit}}$ stable to ±0.5h across that range (matches Methods §Numerical simulation claim).
4. **Commit + tag**: make the v2 commits (methods refinements, repo hygiene, restored files), then `git tag -a v2-submission-2026-05-XX -m "Science Advances aeh5879 v2 submission"` and push tag.
5. **Zenodo fresh deposit**: trigger a new Zenodo deposit at the v2-submission tag via GitHub-Zenodo integration (or manual upload). Capture the new code DOI and dataset DOI.
6. **Fill placeholders**: substitute `[V2_SUBMISSION_SHA]`, `[V2_ZENODO_CODE_DOI]`, `[V2_ZENODO_DATA_DOI]` throughout `main_v2.tex`, `supplementary_v2.tex`, and `cover_letter_v2.tex`. Recompile and verify the values render correctly in the PDF output.
7. **Final reconciliation**: rebuild `numerical_claims_v2.csv` with `commit_sha` column updated to the v2-submission SHA for every row that previously cited `d047d2d` (the multiscale numbers reproduce identically because the model code is unchanged from d047d2d to the v2 SHA).
8. **Verify**: AC reads PDF output of all three documents end to end; spot-checks every number in the manuscript and supplement against the corresponding row in `numerical_claims_v2.csv`.
9. **Submit**: send to Erin with corrected figure files attached, along with a note pointing to the numerical-claims CSV in supplementary materials.

## What's in v2 that wasn't in v1

### Substantive corrections (per disclosure email)
- Replaced phenomenological logistic model with multiscale within-host implementation (Phase 1 tau-leap + Phase 2 ODE + eclipse-phase fixed-delay buffer)
- $t_{\text{crit}}$ values: 16–28h / 68–76h (v1, wrong — these were `seeding_midpoint`) → **34.5h / 60.5h** (v2, verified)
- Compression ratio: ≈3× (v1) → **≈1.75×** (v2, concordant with NHP 1.5–2× empirical range)
- Population bound: 6.8% (v1, wrong — used wrong $t_{\text{crit}}$ and wrong F_access) → three-layer evidence stack (v2): cascade point estimates 0.003–0.7%, outbreak probability 92.7% 10-year, analytical envelope ~11%

### Reframed sections
- **Abstract**: cascade-led γ framing with all three evidence layers
- **Population-level efficacy** (formerly "Population-level efficacy bound"): three-layer evidence stack with explicit cross-citation to companion paper
- **NHP concordance**: honest framing — compression-ratio concordant; absolute timepoints partial; explicit acknowledgment of model's domain limit

### New material
- **Methods §Software environment and reproducibility** — Python version, pinned packages, deterministic seeds, byte-identical reproduction claim
- **Methods §Numerical simulation** — N=500/cell, sensitivity range, multiscale phase descriptions
- **Methods §Sex as a biological variable** — SABV statement
- **Mathematical framework §Assumption: transmission-competent exposure** — U=U scoping (option b)
- **Implications** — U=U sentence (option c)
- **City-stratified analysis §6.1** — composite-indicator-not-proxy text (Tchetgen Tchetgen 2020 + Pearl 2009)
- **Population-level efficacy §Limitation: F_access denominator** — U=U denominator scoping
- **Supplement §S7** — Multiscale within-host model implementation (full algorithmic detail)
- **Supplement §S8** — NHP concordance per-timepoint detail with concordant/non-concordant flags
- **Supplement §S6 worked example** — corrected with verified t_crit and bound math

### Visual changes
- **Figure 1 Panel A**: vertical shaded band at $\log_{10}(VL) < \log_{10}(200)$ marking U=U region (per PARTNER2 / modern operational definition); model not applicable, transmission risk ≈ 0
- **Figure 2 Panel A**: same U=U band
- **Figure 3**: unchanged (AIDSVu data not subsettable by mode of acquisition)

### New citations
- Rodger 2019 (PARTNER2)
- Bavinton 2018 (Opposites Attract)
- Tchetgen Tchetgen 2020/2024 (proximal causal inference)
- Pearl 2009 (Causality)
- Demidont companion analysis (in preparation, cited for cascade and outbreak evidence)

### Title change
- v1 LaTeX: "Route-Specific **and Population-Level** Intervention Limits"
- v2 LaTeX: "Route-Specific Intervention Limits" (matches the submitted v1 docx, which never had "and Population-Level")

## What's NOT changed in v2

- Theorem statements (Theorem S3.1 Finite Prevention Window, Lemma S5.2 Inoculum-Monotone Hitting Times, Corollary S5.1 Route Compression, Proposition S6.1 Population Bound) — all are correct as mathematical statements; no changes.
- Mathematical framework (three-state Markov, master equation, hazard decomposition) — unchanged.
- Equation [3] (city-level structural delay function) — unchanged; only the methodological note about it being a composite indicator (not proxy) is added.
- Equation [4] (population bound formula) — unchanged.
- Figure 3 — unchanged.
- Author name, affiliation (Philadelphia), ORCID, email — unchanged.
- Competing-interests disclosure (Gilead employment / divestiture) — unchanged.

## Notes for Erin's editorial workflow

- Three documents to upload: `main_v2.pdf`, `supplementary_v2.pdf`, `cover_letter_v2.pdf` (compile from the .tex files).
- Figure source files (high-resolution PNGs/PDFs) are attached separately to address the original combined-PDF rendering issue.
- The submission ID `aeh5879` is preserved from the original Science transfer.
- The substantive correction was self-disclosed in prior correspondence; the v2 package implements the corrections with full reproducibility provenance.

## Numerical-claims CSV (companion artifact)

`numerical_claims_v2.csv` (24 rows) consolidates every numerical claim in the v2 manuscript and supplement against its source file, source locator, commit SHA, and verification method. Categories:
- 4 rows: $t_{\text{crit}}$ values (parenteral & mucosal at CV=0.0 and CV=0.3)
- 1 row: compression ratio
- 2 rows: extinction probabilities
- 2 rows: F_access at canonical and JID-consistent medians
- 2 rows: bound under canonical and JID-consistent parameterizations
- 5 rows: NHP concordance per timepoint
- 4 rows: cascade-based point estimates (current policy, decriminalization, MSM comparison)
- 4 rows: outbreak forecasting (national 5/10-yr, regional, median time)

Every row with `commit_sha` value `d047d2d` will be updated to `[V2_SUBMISSION_SHA]` after the v2 tag is created (the multiscale model code is unchanged between d047d2d and the v2 tag, so the numbers reproduce identically; only the cited SHA changes).
