"""
Corner-2 First-Pass Run — BTG v2.2 Barrier Weights (PRE-WISE)
==============================================================

First-pass run of the PEP finite-window layer on literature-grounded
(BTG v2.2.0) barrier weights. Results are labeled PRE-WISE.

WISE-fitted weights will supersede once WISE data are available.
See docs/chen_update_2026-05-31.md in nyx-wise-methods for context.

DO NOT export kinetic constants or PEP numbers to any Chen-visible
or public-facing location. Methods paper integrates by reference;
numbers held until FPW manuscript publishes.
"""

from __future__ import annotations

import json
import sys
import random
import numpy as np
from datetime import datetime
from pathlib import Path

# Path resolution
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from corner2_pep.substrate_fpw import FPW_SUBSTRATE, PEP_ACCESS_WEIGHTS
from corner2_pep.pep_access import BarrierProfile, p_access
from corner2_pep.pep_layer import net_protection, pep_efficacy_at_delay

# ---------------------------------------------------------------------------
# RNG SEED (reproducibility)
# ---------------------------------------------------------------------------
SEED = 20260531
rng = np.random.default_rng(SEED)
random.seed(SEED)

# ---------------------------------------------------------------------------
# BTG v2.2 POPULATION BASELINES (from lai_prep_config.json v2.2.0)
# Loaded by reference — NOT imported from public repo at runtime
# ---------------------------------------------------------------------------
BTG_V22_POPULATIONS = {
    "MSM":               {"baseline_attrition": 0.45, "typical_barriers": ["TRANSPORTATION", "PRIVACY_CONCERNS", "HEALTHCARE_DISCRIMINATION"]},
    "CISGENDER_WOMEN":   {"baseline_attrition": 0.55, "typical_barriers": ["TRANSPORTATION", "MEDICAL_MISTRUST", "INSURANCE_DELAYS"]},
    "TRANSGENDER_WOMEN": {"baseline_attrition": 0.50, "typical_barriers": ["HEALTHCARE_DISCRIMINATION", "MEDICAL_MISTRUST", "PRIVACY_CONCERNS"]},
    "ADOLESCENT":        {"baseline_attrition": 0.65, "typical_barriers": ["TRANSPORTATION", "PRIVACY_CONCERNS", "SCHEDULING_CONFLICTS"]},
    "PWID":              {"baseline_attrition": 0.75, "typical_barriers": ["LEGAL_CONCERNS", "LACK_IDENTIFICATION", "MEDICAL_MISTRUST", "HEALTHCARE_DISCRIMINATION"]},
    "PREGNANT_LACTATING":{"baseline_attrition": 0.55, "typical_barriers": ["TRANSPORTATION", "SCHEDULING_CONFLICTS", "INSURANCE_DELAYS"]},
    "GENERAL":           {"baseline_attrition": 0.47, "typical_barriers": ["TRANSPORTATION", "INSURANCE_DELAYS"]},
    "SEX_WORKER":        {"baseline_attrition": 0.55, "typical_barriers": ["LEGAL_CONCERNS", "HEALTHCARE_DISCRIMINATION", "PRIVACY_CONCERNS", "MEDICAL_MISTRUST"]},
}

# PEP delay scenarios (hours from exposure)
PEP_DELAY_HOURS = [2, 6, 12, 24, 36, 48, 72]


def run_population_sweep() -> list:
    """Sweep all 8 populations × delay scenarios."""
    results = []
    for pop_key, pop_data in BTG_V22_POPULATIONS.items():
        ba = pop_data["baseline_attrition"]
        prep_success = 1.0 - ba
        barriers = pop_data["typical_barriers"]

        for delay_h in PEP_DELAY_HOURS:
            profile = BarrierProfile(
                barrier_keys=barriers,
                averted_by_prep=prep_success,
                population_key=pop_key,
            )
            result = net_protection(
                prep_adjusted_success=prep_success,
                profile=profile,
                hours_to_pep=delay_h,
            )
            results.append({
                "population": pop_key,
                "baseline_attrition": ba,
                "prep_adjusted_success": prep_success,
                "hours_to_pep": delay_h,
                "net_protection": result["net_protection"],
                "pep_gain_over_prep": result["net_protection_gain_over_prep_alone"],
                "pep_deficit": result["pep_deficit"],
                "pep_deficit_pct": result["pep_deficit_pct_of_potential"],
                "eclipse_exceeded": result["eclipse_exceeded"],
                "p_access_raw": result["p_access_raw"],
                "p_access_net": result["p_access_net"],
                "pep_efficacy_at_delay": result["pep_efficacy_at_delay"],
                "active_barriers": result["active_barriers"],
                "label": "PRE-WISE",
            })
    return results


def compute_summary(results: list) -> dict:
    """Aggregate summary statistics across populations and delays."""
    by_pop = {}
    by_delay = {}

    for r in results:
        pop = r["population"]
        delay = r["hours_to_pep"]

        if pop not in by_pop:
            by_pop[pop] = {"net_protection": [], "pep_gain": [], "pep_deficit": []}
        by_pop[pop]["net_protection"].append(r["net_protection"])
        by_pop[pop]["pep_gain"].append(r["pep_gain_over_prep"])
        by_pop[pop]["pep_deficit"].append(r["pep_deficit"])

        if delay not in by_delay:
            by_delay[delay] = {"net_protection": [], "pep_gain": [], "pep_deficit": []}
        by_delay[delay]["net_protection"].append(r["net_protection"])
        by_delay[delay]["pep_gain"].append(r["pep_gain_over_prep"])
        by_delay[delay]["pep_deficit"].append(r["pep_deficit"])

    summary_by_pop = {
        pop: {
            "mean_net_protection": float(np.mean(v["net_protection"])),
            "mean_pep_gain": float(np.mean(v["pep_gain"])),
            "mean_pep_deficit": float(np.mean(v["pep_deficit"])),
        }
        for pop, v in by_pop.items()
    }

    summary_by_delay = {
        str(delay): {
            "mean_net_protection": float(np.mean(v["net_protection"])),
            "mean_pep_gain": float(np.mean(v["pep_gain"])),
            "mean_pep_deficit": float(np.mean(v["pep_deficit"])),
        }
        for delay, v in by_delay.items()
    }

    return {"by_population": summary_by_pop, "by_delay_hours": summary_by_delay}


def write_run_metadata(results: list, output_dir: Path) -> dict:
    """Write run metadata JSON."""
    import subprocess
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(_ROOT.parent),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        commit = "unknown"

    metadata = {
        "run_label": "PRE-WISE",
        "run_date": datetime.now().strftime("%Y-%m-%d"),
        "run_timestamp": datetime.now().isoformat(),
        "rng_seed": SEED,
        "commit_hash": commit,
        "config_source": "BTG v2.2.0 (lai_prep_config.json, feature/v2.2.0)",
        "populations_n": len(BTG_V22_POPULATIONS),
        "pep_delay_scenarios_hours": PEP_DELAY_HOURS,
        "substrate": {
            "tau_eclipse_hours": FPW_SUBSTRATE.tau_eclipse_hours,
            "eclipse_boundary_hours": FPW_SUBSTRATE.eclipse_boundary_hours,
            "reference_log10_vl": FPW_SUBSTRATE.reference_log10_vl,
            "reference_integration_complete_hours": FPW_SUBSTRATE.reference_integration_complete_hours,
            "pep_efficacy_peak": FPW_SUBSTRATE.pep_efficacy_peak,
        },
        "barrier_weights_source": "BTG v2.2.0 PRE-WISE (WISE-fitted will supersede)",
        "access_barriers": list(PEP_ACCESS_WEIGHTS.ACCESS_BARRIER_KEYS),
        "note": (
            "DO NOT export to Chen-visible or public-facing locations. "
            "Methods paper integrates by reference. "
            "Numbers held until FPW manuscript publishes."
        ),
        "n_results": len(results),
    }

    meta_path = output_dir / "run_metadata_prepwise_v1.json"
    meta_path.write_text(json.dumps(metadata, indent=2))
    print(f"  Metadata: {meta_path}")
    return metadata


def main():
    print("=" * 70)
    print("CORNER-2: PEP FINITE-WINDOW LAYER — FIRST-PASS (PRE-WISE)")
    print("=" * 70)
    print(f"  Seed:       {SEED}")
    print(f"  Substrate:  tau_eclipse={FPW_SUBSTRATE.tau_eclipse_hours}h  "
          f"eclipse_boundary={FPW_SUBSTRATE.eclipse_boundary_hours}h")
    print(f"  Barrier weights: BTG v2.2.0 PRE-WISE")
    print()

    # Output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = _ROOT.parent / "results" / f"corner2_prepwise_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run sweep
    print("Running population × delay sweep...")
    results = run_population_sweep()
    print(f"  {len(results)} scenarios computed")

    # Summary
    summary = compute_summary(results)
    print()
    print("Summary by population (mean across delays):")
    for pop, s in summary["by_population"].items():
        print(f"  {pop:<25}  net_prot={s['mean_net_protection']:.3f}  "
              f"pep_gain={s['mean_pep_gain']:.3f}  "
              f"deficit={s['mean_pep_deficit']:.3f}")

    print()
    print("Summary by PEP delay (mean across populations):")
    for delay, s in summary["by_delay_hours"].items():
        print(f"  {delay:>3}h  net_prot={s['mean_net_protection']:.3f}  "
              f"gain={s['mean_pep_gain']:.3f}  "
              f"deficit={s['mean_pep_deficit']:.3f}")

    # Save results
    results_path = output_dir / "results_prepwise_v1.json"
    results_path.write_text(json.dumps({"results": results, "summary": summary}, indent=2))
    print(f"\n  Results: {results_path}")

    # Metadata
    write_run_metadata(results, output_dir)

    # Figure 2
    print("\nGenerating Figure 2...")
    try:
        from corner2_pep.make_figure2 import make_figure2
        make_figure2(results, output_dir)
        print(f"  Figure 2 saved to {output_dir}")
    except Exception as e:
        print(f"  Figure 2 failed: {e}")
        import traceback; traceback.print_exc()

    print()
    print("=" * 70)
    print(f"DONE — PRE-WISE run complete. Output: {output_dir}")
    print("REMINDER: Do not export PEP numbers to Chen-visible or public locations.")
    print("=" * 70)

    return results, summary


if __name__ == "__main__":
    main()
