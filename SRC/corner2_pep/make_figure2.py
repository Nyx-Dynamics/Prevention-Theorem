"""
Figure 2 — One Barrier Model, Two Layers
=========================================

Visualization for Corner-2: shows that the same BTG structural-barrier
profile simultaneously determines:
  - PrEP bridge-period success (Layer 1)
  - PEP in-window access probability (Layer 2)

And their additive net protection + deficit.

Panel layout:
  A. Net protection by population and PEP delay (heatmap)
  B. PEP deficit by population (PWID vs others at 24h and 48h)
  C. One-barrier-model diagram: PrEP + PEP layers fed by same barrier profile

DO NOT export this figure to Chen-visible or public-facing locations.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from typing import List, Dict

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from corner2_pep.substrate_fpw import FPW_SUBSTRATE


POPULATIONS_ORDER = [
    "MSM", "GENERAL", "TRANSGENDER_WOMEN", "CISGENDER_WOMEN",
    "PREGNANT_LACTATING", "SEX_WORKER", "ADOLESCENT", "PWID",
]

DELAYS_PANEL_A = [2, 6, 12, 24, 48, 72]


def make_figure2(results: List[Dict], output_dir: Path) -> None:
    """
    Generate and save Figure 2 (PNG + PDF).

    Parameters
    ----------
    results : list of dicts from run_corner2_prepwise.run_population_sweep()
    output_dir : Path — save location (NOT a public or Chen-visible path)
    """
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.40)

    ax_a = fig.add_subplot(gs[0, :2])   # heatmap — net protection
    ax_b = fig.add_subplot(gs[0, 2])    # deficit bar chart
    ax_c = fig.add_subplot(gs[1, :])    # conceptual two-layer diagram

    _panel_a_heatmap(ax_a, results)
    _panel_b_deficit(ax_b, results)
    _panel_c_diagram(ax_c)

    fig.suptitle(
        "Figure 2 — One Barrier Model, Two Layers (PRE-WISE)\n"
        "BTG v2.2.0 barrier weights; eclipse boundary enforced",
        fontsize=11, y=1.01,
    )

    for ext in ("png", "pdf"):
        out = output_dir / f"Figure_2_PEP_Layer_Corner2_PREPWISE.{ext}"
        fig.savefig(out, dpi=150 if ext == "png" else None,
                    bbox_inches="tight")

    plt.close(fig)


def _panel_a_heatmap(ax, results):
    """Panel A: Net protection heatmap (population × PEP delay)."""
    grid = np.full((len(POPULATIONS_ORDER), len(DELAYS_PANEL_A)), np.nan)

    for r in results:
        pop = r["population"]
        delay = r["hours_to_pep"]
        if pop in POPULATIONS_ORDER and delay in DELAYS_PANEL_A:
            ri = POPULATIONS_ORDER.index(pop)
            ci = DELAYS_PANEL_A.index(delay)
            grid[ri, ci] = r["net_protection"]

    im = ax.imshow(grid, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    plt.colorbar(im, ax=ax, label="Net Protection (PrEP + PEP)")

    ax.set_xticks(range(len(DELAYS_PANEL_A)))
    ax.set_xticklabels([f"{d}h" for d in DELAYS_PANEL_A], fontsize=8)
    ax.set_yticks(range(len(POPULATIONS_ORDER)))
    ax.set_yticklabels([p.replace("_", " ").title() for p in POPULATIONS_ORDER], fontsize=8)

    for ri in range(len(POPULATIONS_ORDER)):
        for ci in range(len(DELAYS_PANEL_A)):
            if not np.isnan(grid[ri, ci]):
                ax.text(ci, ri, f"{grid[ri, ci]:.2f}", ha="center", va="center",
                        fontsize=7, color="black")

    ax.set_title("A. Net Protection (PrEP + PEP Layer)\nby Population × PEP Delay",
                 fontsize=9)
    ax.set_xlabel("Hours from Exposure to PEP Initiation", fontsize=8)

    # Mark eclipse boundary
    eclipse_h = FPW_SUBSTRATE.eclipse_boundary_hours
    matching = [i for i, d in enumerate(DELAYS_PANEL_A) if d >= eclipse_h]
    if matching:
        ax.axvline(x=matching[0] - 0.5, color="red", linestyle="--", linewidth=1.5,
                   label=f"Eclipse boundary ({eclipse_h:.0f}h)")
        ax.legend(fontsize=7, loc="lower right")


def _panel_b_deficit(ax, results):
    """Panel B: PEP deficit at 24h and 48h for key populations."""
    focus_delays = [24, 48]
    focus_pops = ["PWID", "SEX_WORKER", "ADOLESCENT", "MSM", "GENERAL"]

    x = np.arange(len(focus_pops))
    width = 0.35
    colors = ["#d62728", "#ff7f0e"]

    for i, delay in enumerate(focus_delays):
        deficits = []
        for pop in focus_pops:
            match = [r for r in results if r["population"] == pop and r["hours_to_pep"] == delay]
            deficits.append(match[0]["pep_deficit"] if match else 0.0)
        ax.bar(x + (i - 0.5) * width, deficits, width, label=f"{delay}h delay",
               color=colors[i], alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n").replace("SEX\n", "Sex\n").title()
                        for p in focus_pops], fontsize=7)
    ax.set_ylabel("PEP Deficit\n(Potential − Realized)", fontsize=8)
    ax.set_title("B. PEP Access Deficit\nby Population (24h vs 48h)", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.axhline(0, color="black", linewidth=0.5)


def _panel_c_diagram(ax):
    """Panel C: Conceptual two-layer diagram."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")

    # Barrier profile box (center)
    ax.add_patch(plt.Rectangle((3.5, 1.5), 3, 1.2, fill=True,
                                facecolor="#ffffcc", edgecolor="#333", linewidth=1.5))
    ax.text(5, 2.15, "Structural\nBarrier Profile\n(BTG v2.2.0)", ha="center", va="center",
            fontsize=9, fontweight="bold")
    ax.text(5, 1.55, "TRANSPORTATION · LEGAL_CONCERNS\nMEDICAL_MISTRUST · SCHEDULING_CONFLICTS\n...",
            ha="center", va="bottom", fontsize=6.5, color="#555")

    # PrEP layer (left)
    ax.add_patch(plt.Rectangle((0.3, 2.4), 2.5, 1.0, fill=True,
                                facecolor="#c6efce", edgecolor="#375623", linewidth=1.5))
    ax.text(1.55, 2.9, "PrEP Layer\nadj. success rate\n(bridge period)", ha="center",
            va="center", fontsize=8.5)

    # PEP layer (right)
    ax.add_patch(plt.Rectangle((7.2, 2.4), 2.5, 1.0, fill=True,
                                facecolor="#dce6f1", edgecolor="#1f497d", linewidth=1.5))
    ax.text(8.45, 2.9, "PEP Layer\nP_access × efficacy(t)\n(in-window recovery)",
            ha="center", va="center", fontsize=8.5)

    # Arrows from barrier profile to both layers
    ax.annotate("", xy=(2.8, 2.9), xytext=(3.5, 2.15),
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))
    ax.annotate("", xy=(7.2, 2.9), xytext=(6.5, 2.15),
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))

    # Net protection box (bottom center)
    ax.add_patch(plt.Rectangle((2.5, 0.2), 5, 1.0, fill=True,
                                facecolor="#f2dcdb", edgecolor="#943634", linewidth=1.5))
    ax.text(5, 0.7, "Net Protection = PrEP success + PEP recovery\n"
            "Deficit = Potential − Realized (barrier-gated)", ha="center",
            va="center", fontsize=8.5)

    # Arrows from layers to net protection
    ax.annotate("", xy=(3.5, 1.2), xytext=(1.55, 2.4),
                arrowprops=dict(arrowstyle="->", color="#375623", lw=1.5))
    ax.annotate("", xy=(6.5, 1.2), xytext=(8.45, 2.4),
                arrowprops=dict(arrowstyle="->", color="#1f497d", lw=1.5))

    # Eclipse boundary note
    ax.text(5, 0.02, f"Eclipse boundary enforced at {FPW_SUBSTRATE.eclipse_boundary_hours:.0f}h — "
            "PEP efficacy = 0 beyond this point",
            ha="center", va="bottom", fontsize=7.5, color="red", style="italic")

    ax.set_title("C. One Barrier Model → Two Prevention Layers (Conceptual)", fontsize=9)
