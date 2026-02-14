#!/usr/bin/env python3
"""
Reproducibility Suite — Prevention Theorem
===========================================

Reproduces all results and figures from:
  "Finite Prevention Windows Under Irreversible Infection Establishment:
   A Mathematical Framework for HIV Post-Exposure Prophylaxis Timing"

Under review: Journal of Infectious Diseases (JID-84776)
Preprint: DOI 10.20944/preprints202601.1090.v1

Usage:
    python reproduce_all.py

Author: AC Demidont, DO / Nyx Dynamics LLC
"""

import os
import sys
import time
import logging
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "SRC"
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"

sys.path.insert(0, str(SRC))
FIG_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)

LOG_PATH = DATA_DIR / "reproducibility_log.txt"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_PATH, mode='w'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger(__name__)


class Validator:
    def __init__(self):
        self.checks = []

    def check(self, name, condition, detail=""):
        status = "PASS" if condition else "FAIL"
        self.checks.append((name, status, detail))
        sym = "\u2713" if condition else "\u2717"
        log.info(f"  {sym} {name}: {status}  {detail}")
        return condition

    def summary(self):
        passed = sum(1 for _, s, _ in self.checks if s == "PASS")
        return passed, len(self.checks)


def step1_theorem_figures():
    log.info("=" * 70)
    log.info("STEP 1: Prevention Theorem Figures (Fig 1 + Fig 2)")
    log.info("=" * 70)
    os.chdir(str(SRC))
    from prevention_theorem_figures import plot_figure_1, plot_figure_2
    plot_figure_1()
    plot_figure_2()
    os.chdir(str(ROOT))


def step2_route_compression(v):
    log.info("")
    log.info("=" * 70)
    log.info("STEP 2: Route-Specific Compression (JID Manuscript Figures)")
    log.info("=" * 70)
    os.chdir(str(SRC))

    import numpy as np
    from scipy.stats import lognorm
    from route_compression import (
        generate_hitting_times, compute_epep, find_tcrit,
        figure1, figure2
    )

    N, eta = 10000, 0.05
    log.info(f"Generating {N:,} hitting time samples (seed=42)...")
    T_seed_m, T_int_m, T_seed_p, T_int_p = generate_hitting_times(N=N, seed=42)

    log.info(f"  Mucosal:     T_int median = {np.median(T_int_m):.1f}h")
    log.info(f"  Parenteral:  T_int median = {np.median(T_int_p):.1f}h")

    t_arr = np.linspace(0, 200, 1000)
    E_m, _, _ = compute_epep(t_arr, T_seed_m, T_int_m)
    E_p, _, _ = compute_epep(t_arr, T_seed_p, T_int_p)

    tcrit_m = find_tcrit(E_m, t_arr, eta)
    tcrit_p = find_tcrit(E_p, t_arr, eta)
    log.info(f"  t_crit mucosal   = {tcrit_m:.1f}h")
    log.info(f"  t_crit parenteral = {tcrit_p:.1f}h")
    log.info(f"  Compression ratio = {tcrit_m / tcrit_p:.1f}x")

    sigma_ln = np.log(2.0)
    F = lognorm.cdf(tcrit_p, s=sigma_ln, scale=72.0)
    pop_bound = F * 0.95 + (1 - F) * 0.05
    log.info(f"  Population PEP bound = {pop_bound*100:.1f}%")

    log.info("\nValidating manuscript claims...")
    v.check("Mucosal window ~72h", 50 < tcrit_m < 95, f"{tcrit_m:.1f}h")
    v.check("Parenteral window ~12-24h", 8 < tcrit_p < 35, f"{tcrit_p:.1f}h")
    v.check("Compression >2x", tcrit_m / tcrit_p > 2.0, f"{tcrit_m/tcrit_p:.1f}x")
    v.check("E_PEP monotonic (mucosal)", all(E_m[i] >= E_m[i+1] - 1e-10 for i in range(len(E_m)-1)))
    v.check("E_PEP monotonic (parenteral)", all(E_p[i] >= E_p[i+1] - 1e-10 for i in range(len(E_p)-1)))
    v.check("Population bound <10%", pop_bound < 0.10, f"{pop_bound*100:.1f}%")

    pfx = str(FIG_DIR) + "/"
    figure1(E_m, E_p, t_arr, tcrit_m, tcrit_p, save_prefix=pfx + "Fig1_RouteCompression")
    figure2(E_p, t_arr, tcrit_p, save_prefix=pfx + "Fig2_DistributionOverlap")
    os.chdir(str(ROOT))


def step3_pep_mucosal(v):
    log.info("")
    log.info("=" * 70)
    log.info("STEP 3: PEP Mucosal / Reservoir Analysis")
    log.info("=" * 70)
    os.chdir(str(SRC))
    from PEP_mucosal import plot_pep_efficacy_curve, plot_pep_timing_vs_reservoir

    fig1, results = plot_pep_efficacy_curve(save_path=str(FIG_DIR / "pep_efficacy_curve.png"))
    fig2 = plot_pep_timing_vs_reservoir(save_path=str(FIG_DIR / "pep_timing_impact.png"))
    v.check("PEP efficacy figures generated", fig1 is not None and fig2 is not None)
    os.chdir(str(ROOT))


def step4_middle_ground(v):
    log.info("")
    log.info("=" * 70)
    log.info("STEP 4: Continuous Exposure Spectrum Model")
    log.info("=" * 70)
    os.chdir(str(SRC))
    import numpy as np
    from middle_ground import integration_time_params, pep_efficacy_curve
    from middle_ground import plot_spectrum, plot_continuous_window, plot_access_overlap

    V0_values = [1, 10, 100, 1000]
    results = {}
    for V0 in V0_values:
        mu, sigma = integration_time_params(V0)
        hours, efficacy = pep_efficacy_curve(V0)
        results[V0] = {"mu": mu, "sigma": sigma, "hours": hours, "efficacy": efficacy}
        log.info(f"  V0={V0:5d}: median T_int = {np.exp(mu):.1f}h")

    plot_spectrum(results, save_path=str(FIG_DIR / "middle_ground_spectrum.png"))
    plot_continuous_window(save_path=str(FIG_DIR / "middle_ground_continuous.png"))
    plot_access_overlap(save_path=str(FIG_DIR / "middle_ground_overlap.png"))

    mu_1, _ = integration_time_params(1)
    mu_1000, _ = integration_time_params(1000)
    v.check("V0=1 median ~72h", 50 < np.exp(mu_1) < 95, f"{np.exp(mu_1):.1f}h")
    v.check("V0=1000 median ~20h", 10 < np.exp(mu_1000) < 30, f"{np.exp(mu_1000):.1f}h")
    os.chdir(str(ROOT))


def main():
    start = time.time()
    log.info("=" * 70)
    log.info("REPRODUCIBILITY SUITE: PREVENTION THEOREM")
    log.info("Finite Prevention Windows Under Irreversible Infection Establishment")
    log.info(f"Timestamp: {datetime.now().isoformat()}")
    log.info(f"Root:    {ROOT}")
    log.info(f"Figures: {FIG_DIR}")
    log.info(f"Log:     {LOG_PATH}")
    log.info("=" * 70)

    v = Validator()
    failures = []

    steps = [
        ("Prevention Theorem Figures", lambda: step1_theorem_figures()),
        ("Route Compression (JID)", lambda: step2_route_compression(v)),
        ("PEP Mucosal Analysis", lambda: step3_pep_mucosal(v)),
        ("Middle Ground Model", lambda: step4_middle_ground(v)),
    ]

    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            log.error(f"STEP FAILED: {name} — {e}")
            failures.append((name, str(e)))

    elapsed = time.time() - start
    passed, total = v.summary()

    log.info("")
    log.info("=" * 70)
    log.info("REPRODUCIBILITY SUMMARY")
    log.info("=" * 70)
    log.info(f"Validation: {passed}/{total} checks passed")
    log.info(f"Failures:   {len(failures)} step(s) failed")
    log.info(f"Elapsed:    {elapsed:.1f}s")

    if failures:
        for name, err in failures:
            log.info(f"  \u2717 {name}: {err}")

    fig_files = sorted(FIG_DIR.glob("*"))
    if fig_files:
        log.info(f"\nGenerated files ({len(fig_files)}):")
        for f in fig_files:
            log.info(f"  {f.name:45s} ({f.stat().st_size/1024:.0f} KB)")

    if passed == total and not failures:
        log.info("\n\u2605 ALL CHECKS PASSED")
        sys.exit(0)
    else:
        log.info("\n\u26a0 SOME CHECKS FAILED — review above")
        sys.exit(1)


if __name__ == "__main__":
    main()
