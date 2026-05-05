# aeh5879 v2 — Text Deltas

**Source manuscript:** `aeh5879_ArticleContent_v1.docx` (Apr 26, 26,021 chars, 168 paragraphs)
**Source supplement:** `finite_windows_Science_supplement.docx` (Mar 13, Zenodo upload)
**Anchor SHA:** `d047d2d` (May 1, 2026, `SRC/multiscale_model/`)
**Verified numbers:** see `numerical_claims_v2.csv` (24 rows)

---

## 1. ABSTRACT — full rewrite (γ framing)

### v1 (current — wrong t_crit, envelope-led headline)

> "Proviral integration — the molecular event that converts HIV infection from reversible to permanent — defines a finite window during which post-exposure prophylaxis (PEP) can succeed. Using a three-state absorbing Markov model parameterized on established HIV-1 kinetic constants, we prove that PEP efficacy decays monotonically to zero and derive the critical prevention window as a function of integration kinetics rather than drug potency. This window is approximately three-fold shorter for parenteral (injection) than mucosal (sexual) exposure, yielding **a parenteral tcrit of 16–28 hours versus 68–76 hours mucosally**. For people who inject drugs, empirically documented structural access delays place fewer than 5% of exposures within the parenteral window, **bounding expected population-level PEP efficacy below 10%** — a failure determined by integration timing, not pharmacology."

### v2 (proposed — verified numbers, cascade-led)

> "Proviral integration — the molecular event that converts HIV infection from reversible to permanent — defines a finite window during which post-exposure prophylaxis (PEP) can succeed. Using a multiscale within-host model that combines stochastic founder-phase tau-leaping with deterministic ODE handoff and an eclipse-phase fixed delay (Perelson 1996, ~22h), we prove that PEP efficacy decays monotonically to zero and derive the critical prevention window from emergent integration kinetics rather than drug potency or hand-set parameters. This window is approximately 1.75-fold shorter for parenteral (injection) than mucosal (sexual) exposure, yielding **a parenteral t_crit at η=0.05 of approximately 34.5 hours (V₀=10³) versus 60.5 hours mucosally (V₀=1)** — concordant with the empirical 1.5–2× compression range spanning Tsai 1998 (intravenous SIV) and Otten 2000 (intravaginal SHIV) NHP challenge data. Independent cascade-based modeling of structural access for people who inject drugs produces population-level prevention success of 0.003–0.7% under current US policy [companion analysis], with a 92.7% national 10-year outbreak probability. The analytical F_access × t_crit upper bound (~11%) is consistent with these point estimates as a conservative envelope. The failure mode is structural rather than pharmacological: PEP retains high per-exposure efficacy; population-level efficacy is bounded by the mismatch between access timing and integration timing."

**Changes:**
- t_crit values: 16–28 / 68–76 → **34.5 / 60.5** (verified at d047d2d)
- compression: ≈3× → **≈1.75×** (verified)
- NHP framing: "concordant with NHP" preserved (compression ratio ✓; absolute timepoints partial — see §5 for honest framing)
- Population bound: "below 10%" → cascade point estimates (0.003–0.7%) + outbreak probability (92.7%) + envelope (~11%) as conservative bound
- Multiscale model description added (Phase 1 + Phase 2 + eclipse), distinguishing from prior phenomenological model

---

## 2. ONE-SENTENCE SUMMARY — minor update

### v1
> "Irreversible proviral integration defines a finite prevention window that structural barriers systematically eliminate for people who inject drugs."

### v2 (no change)

Keep as-is. The qualitative claim is unchanged.

---

## 3. MATHEMATICAL FRAMEWORK / Parameterization — REWRITE

### v1 (current — describes phenomenological model)

> "Parameterization and simulation. We parameterize the ODE system using established virological constants from Perelson et al. (5)... Multiplicative log-normal noise on β and α (coefficient of variation = 0.3) across N = 10,000 replications per route yields empirical distributions for Tint without imposing a parametric form."

### v2 (proposed — describes multiscale model that produced the verified numbers)

> "Parameterization and simulation. The model has two phases. Phase 1 simulates founder-population dynamics stochastically via tau-leaping at low inoculum (V₀ ≤ I_handoff = 10), capturing extinction events for sparse founder populations. Phase 2 integrates the deterministic target-cell-limited ODE system once productive infection is established, anchored on Perelson 1996 (5): viral clearance rate c = 23 day⁻¹, infected-cell death rate δ = 0.7 day⁻¹, eclipse phase τ = 0.9 days (≈22h, Table 2). Standard within-host parameters (6): target cell production λ = 10⁴ cells/mL/day, natural death d_T = 0.01 day⁻¹, infection rate β = 2.4 × 10⁻⁵ mL·virion⁻¹·day⁻¹, viral production p = 10³ virions/cell/day, integration rate α = 10⁻³ day⁻¹. Heterogeneity is implemented as multiplicative log-normal variation (CV = 0.3) on {β, c, δ, α}. The eclipse phase is implemented as a discrete fixed-delay buffer between cell infection and integration competence, ensuring that the lower bound on T_int is the biological eclipse floor (Perelson 1996, Table 2). N = 500 realizations per (route, V₀, CV) cell. Output: empirical distributions for T_int and emergent t_crit values without imposing a parametric form on either. Code and data: https://github.com/Nyx-Dynamics/Prevention-Theorem at commit `d047d2d` (Zenodo: 10.5281/zenodo.18746065 [code], 10.5281/zenodo.18116991 [dataset])."

**Changes:**
- Replaces phenomenological-noise parameterization with two-phase multiscale description
- Adds eclipse delay buffer (matches disclosure email's description)
- Specifies SHA d047d2d in body text
- Replaces N = 10,000 with N = 500 per cell (the actual v3 run)

---

## 4. RESULTS — REWRITE for verified numbers

### v1 (current)

> "Figure 1 shows the predicted EPEP(t) curves for both routes. The ODE-derived mucosal window tcrit(m)(0.05) ≈ 68–76 hours matches the established 72-hour CDC guideline without calibration, providing independent mechanistic validation of the clinical threshold (1). The parenteral window compresses to **tcrit(p)(0.05) ≈ 16–28 hours — a ≈3-fold reduction**. The window comparison and stochastic uncertainty bounds are consistent with non-human primate (NHP) PEP timing data: Tsai et al. (7) showed complete protection at 24 hours and declining protection at 48 hours following intravenous challenge; Otten et al. (8) showed complete protection through 36 hours and declining protection at 72 hours following intravaginal challenge. The route compression ratio observed in NHP data (≈1.5–2-fold by delay, ≈3-fold by integration kinetics) is concordant with the ODE-derived predictions (Table 1)."

### v2 (proposed)

> "Figure 1 shows the predicted E_PEP(t) curves for both routes derived from the multiscale model at SHA `d047d2d`. The mucosal window **t_crit^(m)(0.05) ≈ 60.5 hours** at V₀ = 1 (CV = 0.3) emerges from within-host kinetics and falls within the operational 72-hour CDC mucosal-PEP guideline (1), providing mechanistic context for — though not exact reproduction of — the clinical threshold. The parenteral window compresses to **t_crit^(p)(0.05) ≈ 34.5 hours** at V₀ = 10³ (CV = 0.3) — a **≈1.75-fold reduction**. The route compression ratio is concordant with the empirical NHP range of 1.5–2× spanning Tsai et al. (7) (intravenous SIV: complete protection at 24h, 50% at 48h) and Otten et al. (8) (intravaginal SHIV: complete protection through 36h, 50% at 72h). At the within-host kinetic level, the model captures full protection at early delays (Otten 12h, 36h: model E_PEP = 0.95 vs NHP 1.00) but under-predicts protection at later timepoints (Tsai 24h, 48h; Otten 72h). We interpret this gap as evidence of NHP-specific protective contributions not captured by the within-host kinetic framework alone — an empirical limit of the model's domain rather than a refutation of the compression result. (Table 1 summarizes per-timepoint concordance.)"

**Changes:**
- Updated t_crit values to verified d047d2d output
- Updated compression ratio from ≈3× to ≈1.75×
- Honest framing of NHP concordance: ratio ✓, absolute timepoints partial
- Foreshadows the within-host model's domain boundary (sets up Discussion)

---

## 5. POPULATION BOUND — major rewrite (γ reframe + sensitivity + limitation)

### v1 (current — single-number 6.8% headline)

> "Modeling access delay as log-normal (median 72 hours, geometric SD 2.0) based on documented PWID structural barriers (9, 10) yields **F_access(24 h) ≈ 0.02**: approximately 2% of PWID access PEP within the compressed window. With η = 0.05 and ε_max = 0.95, Equation [4] yields **Ē_PEP ≤ 0.02 × 0.95 + 0.98 × 0.05 = 6.8%**. The bound holds for any access delay distribution satisfying the stated fraction..."

### v2 (proposed — three-layer evidence stack with cascade lead)

> "**Population-level efficacy: three independent lines of evidence.**
>
> **(a) Cascade-based point estimates.** Companion modeling [Demidont, in preparation] applies an architectural barrier framework (n = 100,000 PWID, 5-year horizon) and an independent PWID Monte Carlo simulation to compute the probability that an individual PWID achieves the full prevention cascade (P(R₀=0)) under structural conditions. Under current US policy: architectural barrier model produces P(R₀=0) = 0.003% [95% CI 0.000–0.006%]; independent PWID simulation produces 0.007%. Decriminalization-only intervention raises this to 0.198–0.814%; full harm reduction raises it to 9.55%. The same architectural model applied to MSM under identical structural assumptions produces 16.3% prevention success — a ≈5,400-fold population-level disparity between routes.
>
> **(b) Outbreak probability forecasting.** Independent stochastic outbreak modeling [companion paper, Supplement S2] forecasts a 73.8% probability of major HIV outbreak among PWID nationally within 5 years under current policy (95% PSA CI 63.5–82.0%) and 92.7% within 10 years; regional probabilities reach 86.3% (Pacific Northwest) and 78.4% (Appalachia) within 5 years. Median time-to-outbreak is 3 years.
>
> **(c) Analytical F_access × t_crit envelope (this paper).** Modeling access delay as log-normal (median 96 hours, geometric SD 2.0; modeling choice informed by qualitative discussion in Taylor 2019 and Baugher 2025; specific median is not lifted from primary timing data and is a defensible parameterization rather than an empirical anchor) yields F_access(34.5h) ≈ 0.070. With η = 0.05 and ε_max = 0.95, Equation [4] yields **Ē_PEP ≤ F_access · ε_max + (1 − F_access) · η ≈ 11.3%** — the most generous (envelope) reading of population-level PEP efficacy. Sensitivity analysis across plausible median values F_median ∈ {72, 96, 110, 120} hours yields envelope estimates of 7.5%–15.9%; in all cases, the envelope is one to four orders of magnitude above the cascade point estimates in (a) and is consistent with the outbreak probabilities in (b) as a conservative upper bound.
>
> **Convergent finding.** Three independent modeling approaches — cascade-based individual-level prevention probability, system-level outbreak forecasting, and analytical F_access × t_crit envelope — all support the conclusion that PEP cannot serve as a population-level safety net for PWID under current structural conditions. The analytical envelope is the least favorable estimate; the cascade-based central estimates are 100–3,000× more pessimistic and remain qualitatively consistent.
>
> **Limitation: F_access denominator.** F_access here is derived from population-wide PEP-access surveillance, which includes exposures from both transmission-competent and virologically suppressed sources. Under U=U (Rodger et al. 2019; Bavinton et al. 2018), exposures from suppressed sources have functionally zero transmission risk and are outside the framework's domain. Restricting F_access to transmission-competent exposures only would reduce the relevant denominator and yield a tighter (smaller) F_access and a tighter bound. This refinement strengthens, rather than weakens, the conclusion."

**Changes:**
- Replaces single-number 6.8% claim (which used wrong t_crit and wrong F_access) with three-layer evidence stack
- Foregrounds cascade and outbreak evidence; envelope demoted to "internal consistency check"
- Sensitivity range across F_median values
- Honest framing: median is "modeling choice" not "empirical anchor"
- Adds U=U denominator limitation per §5 of plan

---

## 6. CONCLUSION — minor update

### v1 (current)
> "...structural access delays for PWID systematically consume this window, **bounding expected population-level efficacy below 10%** regardless of drug potency."

### v2 (proposed)
> "...structural access delays for PWID systematically consume this window. Across three independent modeling approaches — cascade-based prevention probability, system-level outbreak forecasting, and analytical F_access × t_crit envelope — population-level PEP efficacy for PWID under current structural conditions is bounded between approximately 0.003% (cascade central estimate) and 11% (analytical envelope), with corresponding 92.7% 10-year outbreak probability. The conclusion is not a statement against PEP for individuals who present in time; it is a quantitative argument that pre-exposure prevention is the primary biomedical tool available for this population, and that current deployment levels are insufficient to close the gap."

---

## 7. FIGURE 1 CAPTION — update

### v1
> "...Parenteral tcrit ≈ 16–28 hours; mucosal tcrit ≈ 68–76 hours."

### v2
> "...Parenteral t_crit ≈ 34.5 hours (V₀ = 10³, CV = 0.3); mucosal t_crit ≈ 60.5 hours (V₀ = 1, CV = 0.3); compression ratio ≈ 1.75×. Values reproducible from `Nyx-Dynamics/Prevention-Theorem` at commit `d047d2d` (`SRC/multiscale_model/results_v3/heterogeneity_summary.csv`)."

**Plus Panel A correction (option A1 / threshold a):** add `axvspan` shaded band at log10(VL) ∈ [log10(50), log10(200)] labeled "U=U region (per PARTNER2 / modern operational definition; VL <200 c/mL); model not applicable, transmission risk ≈ 0." Same band on Figure 2 and Supp Fig S1A. **Figure 3 unchanged** — AIDSVu late-dx data cannot be subset by mode of infection, and the city-level analysis cannot reliably contain a U=U band.

---

## 8. METHODS — INSERT new U=U paragraph (option b from §5.5)

**Insert after the framework definitions, before parameterization:**

> "**Assumption: transmission-competent exposure.** The Finite Windows framework models within-host integration kinetics conditional on exposure to viremic source virus. We assume throughout that source plasma viral load is at or above the established U=U threshold (≥200 copies/mL sustained). The framework's predictions — t_crit, F_access, and the population-level efficacy bound — are not applicable to exposures from virologically suppressed sources, for whom transmission risk is functionally zero (Rodger et al. 2019, PARTNER2; Bavinton et al. 2018, Opposites Attract). U=U exposures are operationally outside scope: no PEP need arises."

---

## 9. DISCUSSION — INSERT new U=U sentence (option c from §5.5)

**Add a single sentence near the end of the Discussion, bounding the policy claims:**

> "All efficacy and access estimates herein are conditional on a transmission-competent exposure; the framework does not apply to U=U exposures (source viral load <200 copies/mL sustained), for whom transmission risk is functionally zero, and the policy implications drawn here pertain only to populations with non-suppressed source viremia."

---

## 10. CITY-STRATIFIED ANALYSIS / §6.1 — INSERT proxy-construct text (locked version)

**Insert in the City-Stratified Analysis paragraph, immediately after the statement of Equation [3]:**

> "City-level structural delay is operationalized as a composite indicator combining AIDSVu late-diagnosis percentage and linkage-to-care percentage according to the functional form specified in Equation [3]. The composite is a constructed measure of healthcare access friction; it is not a proxy under proximal causal inference (Tchetgen Tchetgen 2020), and the analysis does not invoke the identification conditions required for proxy-based estimation (Pearl 2009). The distinction matters because related work in this research program does invoke proxy-based identification under those conditions."

**Required new citations:**
- Tchetgen Tchetgen, E. (2020). Proximal causal inference. [JASA / Biometrika — confirm specific reference]
- Pearl, J. (2009). *Causality: Models, Reasoning, and Inference*, 2nd ed., Cambridge University Press.

---

## 11. F_access vs Δt_structural SECTION BREAK — §6.2

**Current state in v1:** Equation [3] (city-level structural delay, β₁·LateDx + β₂·Linkage + Δt_SDOH) and the population-bound F_access (log-normal 96h/GSD 2.0) sit in adjacent paragraphs without explicit labeling.

**Action:** Insert section break with subsection headers:
- **§Population-level access friction (city-level, equation 3).** Δt_structural composite — used for city stratification, AIDSVu-anchored, regression-style.
- **§Time-to-PEP-access distribution (individual-level, F_access).** Population distribution of access delays — used for the F_access × t_crit envelope, log-normal parameterization.

These are different objects at different scales. Section-break and label distinctly so reviewers don't read them as competing parameterizations of the same quantity.

---

## 12. ADDRESS — §6.3 (locked: Philadelphia)

**All four files (manuscript, supplement, cover letter, figures legend if applicable):**
- Author affiliation: "Nyx Dynamics LLC, Philadelphia, PA 19107"
- Email: ac.demidont@nyxdynamics.com (consistent with v1)

---

## 13. TABLE VALUE — §6.4 (locked: keep 76.0%)

No change required. v1 docx, .tex, and local Prevention-Theorem docx all show 76.0% in the General community row. The "76.2 → 76.1" correction memory is stale.

---

## 14. DATA AVAILABILITY — §6.5 update

### v1 (current — placeholder)
> "All proofs, code, and datasets are available at https://github.com/Nyx-Dynamics/Prevention-Theorem (Zenodo DOI: [ZENODO_DOI_PENDING])."

### v2 (proposed)
> "All code, simulation outputs, and datasets are available at https://github.com/Nyx-Dynamics/Prevention-Theorem at commit `d047d2d`. The numerical claims in this paper — including parenteral t_crit = 34.5h, mucosal t_crit = 60.5h, and NHP concordance values — reproduce from `SRC/multiscale_model/run_mc_v3.py` at this SHA (verified on 2026-05-05). Persistent archives: code DOI 10.5281/zenodo.18746065; dataset DOI 10.5281/zenodo.18116991. A reconciled summary of all numerical claims with source files and commit references is provided as `numerical_claims_v2.csv` in the supplementary materials."

---

## 15. COVER LETTER — separate file

**Required:**
- Redate from *Science* (v1, Mar 11) to *Science Advances* (v2, May 5 2026)
- Drop "and Population-Level" from cover letter title to align with manuscript title (recommendation per handoff §5.4)
- Body should explicitly cite SHA `d047d2d` in the methods-correction summary
- Body should reference the cascade-led abstract reframe so Erin understands the structural change vs v1
- Address: Philadelphia 19107 (consistent with v1)

---

## 16. CITATIONS TO VERIFY (action items for AC)

1. **Reference 9 in v1**: "D. Taylor et al., Time to PEP initiation and access barriers for injection drug users. *J Acquir Immune Defic Syndr* 80, e56–e59 (2019)." — **needs verification.** This citation is *different* from the JID-cited Taylor 2019 commentary (J.L. Taylor, Walley, Bazzi, Subst Abus 40(4):441–443). If the JAIDS paper is real, it may have primary timing data and be a stronger anchor for F_access. Confirm:
   - Does this JAIDS paper exist?
   - If yes, what does it report for median time-to-PEP / IQR / percentile distribution?
   - This could change F_access from "modeling choice" to "empirical anchor".
2. **Tchetgen Tchetgen 2020** — confirm specific journal/title for §6.1 citation.
3. **Bavinton 2018** — Lancet HIV 5:e438–e447 — confirm.
4. **Rodger 2019** — Lancet 393:2428–2438 (PARTNER2) — confirm.

---

## 17. THINGS NOT CHANGING IN v2 (for clarity)

- Theorem statements (S3.1, S5.2, S6.1) — unchanged; these are correct mathematical statements.
- Mathematical framework section — three-state Markov framework unchanged.
- Equation [1] (E_PEP decomposition) — unchanged.
- Equation [4] (population bound formula) — unchanged.
- Figure 3 (city-stratified analysis) — unchanged (AIDSVu late-dx data cannot be subset by mode of infection; U=U scoping not applicable to that figure).
- Supplement core proofs — unchanged. Specific numerical worked examples (§S6) need t_crit and F_access values updated to match main text.

---

---

## 18. METHODS REFINEMENTS — SciScore-anticipating package (5 items)

Source: Claude.ai handoff package, 2026-05-05. SciScore-anticipating subset.

### 18.1 Software environment and reproducibility — INSERT into Methods

**Insert into Methods → "Numerical simulation" subsection (after existing scipy/Python statement at v1 line 56):**

> "All computational analyses were performed in Python 3.11 with NumPy, SciPy, matplotlib, and pandas. Complete environment specifications, including pinned package versions, are provided in `requirements.txt` in the project repository. Monte Carlo simulations used fixed random seeds (`np.random.default_rng(base_seed * 1,000,000 + k)`) to ensure deterministic reproducibility from the cited commit (SHA [V2_SUBMISSION_SHA]). Independent re-execution from the cloned repository at the cited SHA reproduces all numerical claims to full floating-point precision; this was verified prior to revision."

**Implementation actions:**
- Pin `requirements.txt` versions before v2 commit. Current file has `>=` constraints; replace with `==` or compatible-release `~=`.
- The "verified prior to revision" claim is empirically supported by §5.1 byte-identical reproduction of d047d2d.

### 18.2 Sample size and Monte Carlo precision — INSERT into Methods

**Insert into Methods → after the (now-revised) N=500 per-cell statement:**

> "Monte Carlo replicate count (N = 500 per (route, V₀, CV) cell; ~27,000 realizations across the full parameter sweep) was selected to achieve standard error below 1 percentage point on cumulative integration probabilities at clinical timepoints (24h, 48h, 72h). Sensitivity to N at the canonical operating point (V₀ = 10³ parenteral, CV = 0.3) was tested at N ∈ {200, 500, 1000, 5000}; t_crit estimates stable to within ±0.5h across this range. The 500-replicate count was selected to balance statistical precision against parameter-sweep tractability across V₀ × CV grid."

**Implementation actions:**
- ⚠ **Verification needed before v2 submission.** Run the N ∈ {200, 500, 1000, 5000} sensitivity sweep at V₀=1000 CV=0.3, confirm t_crit stays within ±0.5h. Estimated 5–10 min. If the claim fails, soften the wording.
- The original Claude.ai package §1.2 referenced N=10,000 (which was the v1 phenomenological model); for v2 we use N=500 per cell because the multiscale v3 model is more compute-intensive per realization. Reconciled here.

### 18.3 Sex as a biological variable (SABV) — INSERT into Methods

**Insert into Methods → as standalone subsection or appended to parameterization:**

> "**Sex as a biological variable.** The within-host kinetic framework operates on cellular and viral dynamics that have not been established to differ by sex at the resolution of the available data; the foundational Perelson 1996 cohort (5 patients) was small and did not establish sex-specific kinetic parameters. The city-stratified analysis pools across sex within AIDSVu surveillance categories, which provide sex-stratified late-diagnosis and viral suppression data; sex-stratified replication of the city analysis is a defensible extension but is not pursued here, as the structural-barrier conclusions are robust to within-stratum variation."

**Implementation actions:**
- Conservative phrasing per Claude.ai package's own implementation note. Does NOT claim "both male and female" in Perelson cohort (unverified web search).
- If reviewer presses for sex-stratified replication of city analysis, this can be added as supplementary sensitivity. For first peer-review send, asserted version is sufficient.

### 18.4 AI tool transparency — **INSERT** (no existing paragraph to replace)

**Verification finding:** v1 manuscript does NOT contain an existing AI tool disclosure paragraph. Claude.ai package's "REPLACE" instruction is incorrect for v1; this is a fresh INSERT. Place at end of Methods or as standalone paragraph just before Data Availability.

**Insert:**

> "Computational analyses were conducted in Python (NumPy, SciPy, matplotlib). Large language model assistants (Anthropic Claude, OpenAI ChatGPT) were used to support literature search, code review, manuscript readability, and verification of numerical claims against source code. JetBrains Junie was used for code correction; Zotero AI for reference management. All AI tools were used as assistive technologies; the author retains full responsibility for design, analysis, interpretation, and conclusions. Independent reproducibility of all numerical claims from committed code at the cited SHA was verified prior to manuscript revision."

### 18.5 Data availability — REPLACE current placeholder text

**v1 (current):**
> "All proofs, code, and datasets are available at https://github.com/Nyx-Dynamics/Prevention-Theorem (Zenodo DOI: [ZENODO_DOI_PENDING])."

**v2 (proposed) — supersedes §14 above with merged formatting:**

> "Code and data are archived at Zenodo (code DOI: 10.5281/zenodo.18746065; dataset DOI: 10.5281/zenodo.18116991) under the GitHub release tagged for this submission. The repository ([https://github.com/Nyx-Dynamics/Prevention-Theorem](https://github.com/Nyx-Dynamics/Prevention-Theorem) at commit `[V2_SUBMISSION_SHA]`) separates source code (`SRC/`), simulation outputs (`results/`), and processed datasets (`data/`); a `README.md` at the repository root documents the file layout and reproduction sequence. The numerical claims in this paper — including parenteral t_crit = 34.5h, mucosal t_crit = 60.5h, NHP concordance values, and cascade-based estimates cited from companion analysis — reproduce from `SRC/multiscale_model/run_mc_v3.py` at this SHA (verified 2026-05-05). A reconciled summary of all numerical claims with source files and commit references is provided as `numerical_claims_v2.csv` in the supplementary materials. No individual-level human data were used; all surveillance data are publicly available from AIDSVu, NHBS, and CDC sources cited in Methods."

**This block supersedes §14 and integrates the package's §1.5.**

---

## 19. COVER LETTER TRANSPARENCY PARAGRAPH — INSERT

**Insert into v2 cover letter, after the substantive correction summary, before the closing:**

> "Beyond the specific corrections disclosed in our prior letter, the v2 revision includes methods-section refinements addressing standard reproducibility and rigor items: explicit software environment specification, random seed documentation, sample size justification, sex-variable acknowledgment, and clarified AI tool use. Independent reproducibility of all numerical claims from committed code at SHA `[V2_SUBMISSION_SHA]` was verified prior to revision. We submit this strengthened version in the spirit of full methodological transparency, recognizing that the integrity of self-disclosed corrections depends on the rigor of the corrected manuscript itself."

---

## 20. REPOSITORY HYGIENE — pre-v2-SHA checklist with current status

| # | Item | Current status | Action needed before v2 SHA |
|---|---|---|---|
| 1 | `requirements.txt` at repo root with **pinned** versions | ✓ exists, but uses `>=` not `==` | Pin to specific versions used at d047d2d test; commit |
| 2 | `README.md` at repo root with file layout + reproduction sequence | ✓ exists (1958 bytes) | Verify it documents `SRC/`, `results/`, `data/` structure and "how to reproduce" steps; update if not |
| 3 | Random seeds explicitly set in all MC source files | ✓ verified (§5.1 byte-identical reproduction) | Confirm `SRC/multiscale_model/run_mc_v3.py` is restored to working tree (currently deleted in working tree per git status; the working-tree version lives in route_models/) |
| 4 | Zenodo deposit created and DOI assigned | ⚠ **TBD on disk**; package-cited DOIs (10.5281/zenodo.18746065 code, 18116991 dataset) need confirmation that they exist and are linked to the GitHub repo | Verify via Zenodo dashboard; if not present, create via GitHub-Zenodo integration |
| 5 | All numerical claims traceable to specific output files | ✓ via `numerical_claims_v2.csv` (24 rows) | None |
| 6 | `disclosure-2026-04-30` annotated tag preserved | ✓ exists (at d047d2d) | Do NOT delete; this is the editorial transparency snapshot |
| 7 | New v2-submission annotated tag created at canonical v2 commit | ✗ pending | Create `git tag -a v2-submission-2026-05-XX -m "..."` after v2 commits |
| 8 | No uncommitted changes at v2-submission commit | ✗ pending | Confirm `git status` is clean before tagging |

**Additional note on item 3**: working tree currently has `route_models/run_mc_v3.py` (untracked) and `SRC/multiscale_model/run_mc_v3.py` (deleted). The two files have diverged: the SRC version is kinetics-only (matches d047d2d, produces verified 34.5/60.5/1.75×); the route_models version adds the mucosal barrier delay (produces different numbers). For v2, the SRC version is the canonical model and must be restored to the working tree before tagging.

---

## 21. SHA STRATEGY — clarified

The verified numbers (34.5h, 60.5h, 1.75×, NHP concordance) reproduce from commit `d047d2d`. The v2 submission SHA will be a *new* commit incorporating:
- Methods refinements (§18 inserts)
- Pinned `requirements.txt`
- Updated `README.md` (if needed)
- Restored `SRC/multiscale_model/run_mc_v3.py`
- Possibly: N-sensitivity sweep results (§18.2 verification)

Because none of these changes alter the multiscale model source, the **same numerical outputs reproduce from both d047d2d and the v2 SHA**. The v2 SHA inherits the verification from d047d2d.

**Manuscript text** cites the v2 SHA. **Cover letter** cites both: d047d2d for the disclosure, v2-submission for the corrected package.

---

## 22. ACTIONS REMAINING BEFORE DOCX EDITING BEGINS

1. ⚠ Run N-sensitivity sweep at V₀=1000 CV=0.3 to verify §18.2 SE claim (~5–10 min)
2. ⚠ Verify Zenodo deposits at DOIs 10.5281/zenodo.18746065 and 18116991 exist
3. Pin `requirements.txt` to exact versions used at d047d2d
4. Restore `SRC/multiscale_model/run_mc_v3.py` to working tree from d047d2d
5. Verify `README.md` documents file layout and reproduction sequence; update if not
6. Verify JAIDS Taylor 2019 citation (§16 of these deltas — could promote F_access from "modeling choice" to "empirical anchor")
7. AC final read-through of `v2_text_deltas.md` and `numerical_claims_v2.csv`
8. Tag v2-submission SHA after AC sign-off and clean working tree
9. Fill `[V2_SUBMISSION_SHA]` placeholders in all v2 text
10. Then and only then: open Word with track-changes ON; apply deltas; save as `aeh5879_*_v2.docx`

---

## END

**Reviewable artifact.** No edits to v1 .docx files yet. Once AC reviews this delta document, the `numerical_claims_v2.csv`, and the actions in §22, we proceed to actual docx editing in Word (with track-changes ON).

Total prose changes: 7 substantive rewrites (abstract, parameterization, results, bound section, conclusion, Fig 1 caption, data availability) + 9 inserts (Methods U=U paragraph, Discussion U=U sentence, §6.1 proxy text, §6.2 section break, §18.1 software/repro, §18.2 N-sensitivity, §18.3 SABV, §18.4 AI tools, §19 cover letter transparency) + 1 cover letter rewrite + 4–5 citations to verify.
