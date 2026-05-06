# Numerical Claims Verification Report — v2 Submission (RE-RUN)

**Date:** 2026-05-06 05:12:02
**v2-submission SHA:** `dfff300`

## Headline

**24/24 claims verified.** AC adopted Option A on 2026-05-06: Table 3 sensitivity rows recomputed at canonical t_crit=34.5h. Manuscript narrative and CSV are now consistent.

## Per-claim summary

| claim_id | expected | computed | match |
|---|---|---|---|
| `t_crit_par_v3_cv03` | 34.5 | 34.5 | ✓ |
| `t_crit_muc_v3_cv03` | 60.5 | 60.5 | ✓ |
| `t_crit_par_v3_cv0` | 32.5 | 32.5 | ✓ |
| `t_crit_muc_v3_cv0` | 57.5 | 57.5 | ✓ |
| `compression_ratio_v3_cv03` | 1.75 | 1.75 | ✓ |
| `p_extinct_muc_v3_cv03` | 68.4 | 68.4 | ✓ |
| `p_extinct_par_v3_cv03` | 0.0 | 0.0 | ✓ |
| `F_access_at_t_crit_par` | 0.0699 | 0.0699 | ✓ |
| `F_access_at_t_crit_par_alt` | 0.144 | 0.144 | ✓ |
| `bound_envelope_par` | 11.3 | 11.3 | ✓ |
| `bound_envelope_par_jid_consistent` | 18.0 | 18.0 | ✓ |
| `tsai_24h_cv03` | 0.50 vs 1.00 | 0.50 vs 1.00 | ✓ |
| `tsai_48h_cv03` | 0.00 vs 0.50 | 0.00 vs 0.50 | ✓ |
| `otten_12h_cv03` | 0.95 vs 1.00 | 0.95 vs 1.00 | ✓ |
| `otten_36h_cv03` | 0.95 vs 1.00 | 0.95 vs 1.00 | ✓ |
| `otten_72h_cv03` | 0.0158 vs 0.50 | 0.0158 vs 0.50 | ✓ |
| `arch_current_policy_R0_zero` | 0.003 | 0.003 | ✓ |
| `pwid_sim_current_policy_R0_zero` | 0.007 | 0.007 | ✓ |
| `pwid_sim_decrim_R0_zero` | 0.814 | 0.814 | ✓ |
| `arch_msm_comparison` | 16.3 | 16.3 | ✓ |
| `outbreak_5yr_national` | 73.8 | 73.8 | ✓ |
| `outbreak_10yr_national` | 92.7 | 92.74 | ✓ |
| `outbreak_median_yrs` | 3.0 | 3.0 | ✓ |
| `outbreak_pnw_5yr` | 86.3 | 86.35 | ✓ |


## Resolution log

- 2026-05-05 (initial run): 22/24 verified, 2 flagged.
  - F_access_at_t_crit_par_alt (CSV: 0.121) and bound_envelope_par_jid_consistent (CSV: 15.9%) reproduced only at t_crit=32h, not the canonical t_crit=34.5h cited in main_v2.tex Table 3 caption.
- 2026-05-06 (AC decision): adopted Option A. Table 3 + supplement S6.1 + CSV updated at canonical t_crit=34.5h. New values: F_access(72h, 34.5h) = 0.144; bound (median 72h) = 18.0%. Sensitivity range now 8.2-18.0%.
- 2026-05-06 (this re-run): 24/24 verified.

## Verdict

✓ **All 24 numerical claims now verified.** v2 submission is numerically clean.
