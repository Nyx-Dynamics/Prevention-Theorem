# v3 Repository Build — Handoff to Claude-in-AIR

**Repo:** github.com/Nyx-Dynamics/Prevention-Theorem
**Current state:** v2-prep-2026-05-05 (commit 37e27ea), v2-submission released, README v2-focused
**Target state:** v3 ready for PLOS Computational Biology submission
**Scope of this handoff:** Add R3 PK-driven framework and S14 pharmacy sensitivity to the repo as the v3 contribution layer, build a v3 reproducibility script, update README for PLOS Comp Bio shop, tag and Zenodo-release.

---

## Current repo state (verified 2026-06-03)

### Existing directory structure

```
Prevention-Theorem/
├── SRC/
│   ├── multiscale_model/        # Phase 1/2 within-host model (v2 canonical)
│   ├── perelson/                # Stochastic PEP efficacy + VL uncertainty
│   └── ...
├── aidsvu datasets/             # AIDSVu county/state surveillance
├── archive_exploratory/         # Stale v1 artifacts
├── core_theorem/                # Mathematical core (theorems, proofs)
├── results/                     # Output directory (city_analysis subdirectory)
├── route_models/                # Route-dependent compression scaffolding
├── stochastic_layers/           # Stochastic envelope analyses
├── v2_revision/                 # v2 submission artifacts (figures, scripts, claims CSV)
├── reproduce_all_v2.py          # v2 master reproduction script
├── audit_check.py               # Reproducibility audit
├── audit_verification_output_2026-04-26.txt
├── Makefile                     # `make all` reproduces v2
├── README.md                    # v2-focused; needs v3 update
├── CITATION.cff                 # Citation metadata
├── LICENSE                      # MIT
└── requirements.txt             # Pinned dependencies
```

### Existing release / tag history

- `v1.3.0` — Science Submission v5 — aeh1546 (Mar 14, 2026)
- `disclosure-2026-04-30` — editorial-transparency snapshot at the time of the voluntary numerical-error disclosure (commit `d047d2d`)
- `v2-prep-2026-05-05` — the SHA cited in the v2 manuscript Methods (commit `37e27ea`)
- `v2-submission-2026-05-05` — polished SA submission state

### What's NOT in the repo yet

- `v3_revision/` directory (does not exist on `main`)
- `v3_revision/r3_pk_pd/` — PK-driven framework code (5 files, ~48 KB, see Section 2 below)
- `v3_revision/data/Table_34_cities_full.csv` — 34-city input data for pharmacy sensitivity
- `v3_revision/numerical_claims_v3.csv` — verified v3 numerical claims
- `reproduce_all_v3.py` — v3 master reproduction script (see Section 2)
- `README.md` updates for PLOS Comp Bio submission framing
- v3 release tag and Zenodo deposit

---

## What needs to happen (in order)

### Step 1 — Add the R3 PK-driven framework code

Create directory `v3_revision/r3_pk_pd/` and commit five Python files:

| File | Purpose | Size |
|------|---------|------|
| `pk_model.py` | One-compartment Bateman PK with intracellular metabolite compartments (TFV-DP, FTC-TP, DTG plasma). Parameter values anchored to published clinical PK literature. | 14.2 KB |
| `pd_model.py` | Hill-function pharmacodynamics on normalized active fraction, multiplicative survival across drug classes. | 8.5 KB |
| `effective_epsilon.py` | Computes ε_drug(t_PEP) as the time-averaged combined suppression over the remaining route-specific window [t_PEP, t_crit]. Computes p_clear(Z) stage-dependent biological clearability. | 6.9 KB |
| `pharmacy_sensitivity.py` | S14 sensitivity sweep. Reconstructs per-city structural delay from late_dx_pct and retention; layers Δt_pharmacy ∈ {0, 2, 4, 6, 8, 10, 12} h; recomputes envelope bound and per-city E_PEP. | 9.7 KB |
| `test_regression.py` | Verifies PK-driven framework recovers multiscale baseline (t_crit, envelope) within calibration tolerance at perfect adherence. | 8.8 KB |

All five files have been **path-cleaned for the repo** (no `/mnt/`, `/home/claude/`, or `/Users/` hardcoded). They use `Path(__file__).resolve().parents[2]` to resolve repo root at runtime.

### Step 2 — Add the 34-city input data

`v3_revision/data/Table_34_cities_full.csv` (41 lines, ~6 KB). This is the AIDSVu-derived input that `pharmacy_sensitivity.py` reads. Columns:

```
city, state, late_dx_pct, gamma_enrolled, retention,
Omega_star, deflation_pct, B_IRR, IRR_attenuation_pct
```

Note: this file already exists in the project workspace; AC has the canonical copy. AIR should verify before committing that the workspace version is current.

### Step 3 — Add `reproduce_all_v3.py` at repo root

Mirrors the `reproduce_all_v2.py` pattern with the v3 additions:

1. **Part 1 — Replays v2 main-text figures (1, 2, 3, S1)** by calling the existing v2 scripts. v3 main-text figures are identical to v2 because the PK-driven framework recovers the baseline at perfect adherence; v2 reproduction is the canonical regenerator.
2. **Part 2 — v3-specific analyses:**
   - `r3_pk_pd/test_regression.py` — validates PK-driven against multiscale baseline
   - `r3_pk_pd/pharmacy_sensitivity.py` — generates S14 sweep
3. **Part 3 — Numerical claims verification** — copies `v3_revision/numerical_claims_v3.csv` to the timestamped run directory and prints the eight headline v3 claims.

Output structure: `v3_revision/runs/run_<timestamp>/` containing all figures + CSVs from the run.

### Step 4 — Generate `v3_revision/numerical_claims_v3.csv`

The v3 manuscript references this. It should be the v2 claims CSV (already in repo) **plus** the v3-specific additions. v3 numerical claims that the manuscript depends on, all of which the reproduction pipeline must reproduce:

| Claim | Value | Source |
|------|------|--------|
| Parenteral t_crit at η=0.05 (V0=10³, CV=0.3) | 34.5 h | Multiscale baseline (carries from v2) |
| Mucosal t_crit at η=0.05 (V0=1, CV=0.3) | 60.5 h | Multiscale baseline (carries from v2) |
| Compression ratio (mucosal / parenteral) | ~1.75× | Derived from above |
| Eclipse phase delay | ~22 h | Perelson 1996 anchor |
| F_access × t_crit envelope bound (PWID, lognormal median=96h, GSD=2.0) | ~11.3% | v2 envelope analysis |
| PK-driven t_crit at ρ=1.0 (parenteral) | 34.0 h (Δ = 0.5 h) | r3_pk_pd/test_regression.py |
| PK-driven t_crit at ρ=1.0 (mucosal) | 59.5 h (Δ = 1.0 h) | r3_pk_pd/test_regression.py |
| PK-driven envelope at ρ=1.0 | 11.2% (Δ = 0.1 pp) | r3_pk_pd/test_regression.py |
| Hartford E_PEP collapse at Δt_pharmacy = 8 h | 0.47 → 0.00 | r3_pk_pd/pharmacy_sensitivity.py |
| Envelope across pharmacy sweep | 11.5% → 9.6% | r3_pk_pd/pharmacy_sensitivity.py |
| Mean E_PEP across 34 cities across pharmacy sweep | 96.1% → 83.7% | r3_pk_pd/pharmacy_sensitivity.py |

AIR should generate this CSV by running the test_regression and pharmacy_sensitivity scripts and capturing their outputs. Format should match v2 claims CSV structure.

### Step 5 — Update README.md for v3 / PLOS Comp Bio

Replace the v2-focused README content. Key changes:

- **Title block:** Change SA submission framing to PLOS Comp Bio. New manuscript title: *"The HIV Post-Exposure Prophylaxis Window: A Multiscale Framework Linking Within-Host Integration Kinetics to Population-Scale Structural Access"*
- **Subtitle:** *"PLOS Computational Biology submission"* (drop the SA aeh5879 ID; that submission is closed)
- **Quickstart:** Add `python3 reproduce_all_v3.py` as the canonical command; keep `reproduce_all_v2.py` available as the baseline-reproduction path
- **File inventory:** Add `v3_revision/` row, including `r3_pk_pd/` and the new master script
- **Order of operations:** Document the two-part v3 pipeline (v2 base → R3/S14 layer)
- **Tags section:** Add the planned v3 release tag (see Step 7)
- **Zenodo DOI:** Add a placeholder line for the v3 Zenodo deposit; AIR should NOT mint the DOI until AC confirms the submission package is ready

Suggested edits-only diff approach so the README history shows v2 → v3 transition cleanly rather than a from-scratch rewrite.

### Step 6 — Update `audit_check.py` (or create v3 variant)

The existing `audit_check.py` verifies v2 numerical claims against `v2_revision/numerical_claims_v2.csv`. Two options:

- **Option A (lower-touch):** Create `audit_check_v3.py` that mirrors the existing audit and additionally validates v3 claims (PK-driven recovery, pharmacy sensitivity)
- **Option B (cleaner):** Parameterize `audit_check.py` to accept a `--version v2|v3` argument and read the corresponding claims CSV

Recommend Option A for incremental change; Option B is the architectural-cleanup target for a future PR.

### Step 7 — Tag and Zenodo-release v3

After Steps 1–6 are committed and pushed:

1. Create branch `v3-prep-2026-06-XX` (use the actual prep date)
2. Open PR to main
3. After merge, tag the merged commit as `v3-prep-2026-06-XX` (the commit SHA cited in the v3 manuscript Methods)
4. Trigger Zenodo to mint a new version under DOI `10.5281/zenodo.20044747` family
5. Update the v3 manuscript Methods section with the actual prep tag commit SHA (this replaces the current `37e27ea` reference, which is v2)
6. After PLOS Comp Bio submission, tag a `v3-submission-2026-06-XX` snapshot

**DO NOT push or tag until AC confirms.** The v3 manuscript currently references commit `37e27ea` (v2) for reproducibility — AC may want to keep that reference as the canonical multiscale-baseline commit and use a separate v3 commit for the PK-driven additions, or may want a single new commit. That's a decision AC needs to make explicitly before tagging.

---

## File package delivered to AIR

```
v3_handoff/
├── HANDOFF_v3_repo_build.md            # this document
├── reproduce_all_v3.py                  # master reproduction script
├── r3_pk_pd/
│   ├── pk_model.py                      # path-cleaned
│   ├── pd_model.py                      # path-cleaned
│   ├── effective_epsilon.py             # path-cleaned
│   ├── pharmacy_sensitivity.py          # path-cleaned (uses Path(__file__).resolve().parents[2])
│   └── test_regression.py               # path-cleaned
└── data/
    └── Table_34_cities_full.csv         # 34-city input data
```

Place these in the repo at:
- `r3_pk_pd/*` → `v3_revision/r3_pk_pd/*`
- `data/Table_34_cities_full.csv` → `v3_revision/data/Table_34_cities_full.csv`
- `reproduce_all_v3.py` → repo root

---

## Validation criteria — what AIR should verify before opening the PR

After Steps 1–4, AIR should be able to run from a fresh clone:

```bash
pip install -r requirements.txt
python3 reproduce_all_v3.py
```

and produce, in the timestamped run directory under `v3_revision/runs/`:

- Figure_1_VL_Compression_PEP_Parenteral.{png,pdf}
- Figure_2_Stochastic_Efficacy_VL_Uncertainty.png
- Figure_3_City_Stratified_PEP_Efficacy.png
- Figure_S1_CityComparison_Focus.png
- pharmacy_sensitivity_results.csv
- numerical_claims_v3.csv
- All supporting CSVs (mc_efficacy_curves, mc_realizations, mc_summary, city_vl_profiles, pep_stochastic_*, etc.)

Numerical reproduction (must match to within stated tolerance):

- Parenteral t_crit = 34.5 ± 0.5 h
- Mucosal t_crit = 60.5 ± 0.5 h
- Envelope bound = 11.3% ± 0.5 pp
- PK-driven t_crit (parenteral) at ρ=1.0 = 34.0 ± 0.5 h
- PK-driven envelope at ρ=1.0 = 11.2% ± 0.5 pp
- Hartford E_PEP at Δt_pharmacy = 8 h: collapse to < 0.05
- Envelope monotonically decreasing across pharmacy sweep

If any of these fail, **stop and flag to AC** before proceeding to PR/tag/release. Numerical drift between local PyCharm execution and the manuscript-reported values is the most common cause of reproducibility complaints from reviewers; catching it pre-submission is the point of this verification step.

---

## Things AIR should NOT do without AC's explicit go-ahead

1. **Do not push to main directly.** All changes go through a feature branch and PR. The v3 manuscript is in cover-letter prep; commit history should be clean.
2. **Do not mint the Zenodo DOI** until AC confirms the submission package (manuscript, supplement, cover letter, code) is finalized. DOI minting is irreversible.
3. **Do not delete or restructure the v2_revision/ directory.** v2 reproducibility remains a required path (for the prior-review disclosure trail and for the baseline reproduction in `reproduce_all_v3.py` Part 1).
4. **Do not modify the commit `37e27ea` reference in the v3 manuscript Methods** — that's AC's editorial call.
5. **Do not engage with any non-AC accounts on the repo** (no responding to issues, no replying to comments). The repo is currently public but AC-controlled.

---

## Open questions for AC

1. **Commit strategy for v3:** Single commit on a `v3-prep-YYYY-MM-DD` branch, or split into multiple logical commits (one for r3_pk_pd/, one for pharmacy data, one for reproduce_all_v3.py, one for README)? Multiple commits gives a cleaner git log; single commit is simpler for the Zenodo deposit.
2. **Whether to update the manuscript's cited SHA to a new v3 commit, or keep `37e27ea` and add a "v3 additions" reference.** The current v3 manuscript still cites `37e27ea` as the canonical commit. If AIR creates a new v3 commit and tags it, the manuscript should be updated to cite that new SHA — but only after AC reviews and approves.
3. **Whether to archive v1 references (Lancet HIV, Science, Science Advances) anywhere in the repo metadata.** Current README mentions SA explicitly; v3 README should reflect that the SA submission is closed. Whether to mention the prior submission path in the README at all is a curation choice.
4. **Whether AIR should run a final compile of main + supplement LaTeX inside the repo as a sanity check** before tagging. The current v3-3 / v3-4 TeX files compile in AC's local environment; AIR should verify the same against the repo state.
