"""
envelope_corridor_figure.py — Figure S14 (v4)

Two-panel layout:
  Panel A (left, wider):  Envelope corridor with 34 high-burden US metros
                          overlaid at their Δt_pharm=0 baseline positions.
  Panel B (right, zoom):  Hartford pharmacy displacement trajectory across
                          Δt_pharm ∈ {0, 2, 4, 6, 8, 10, 12} hours, showing
                          the cliff crossing at Δt_pharm=10h.

Reads:
  v3_revision/results/pharmacy_sensitivity_corrected/envelope_corridor.csv
  v3_revision/results/pharmacy_sensitivity_corrected/city_envelope_positions.csv

Produces:
  v3_revision/results/pharmacy_sensitivity_corrected/Figure_S14_envelope_corridor.{png,pdf}

Design notes
------------
- Two panels rather than one let each story breathe (corridor geometry on
  the left, displacement mechanics on the right).
- Legends are panel-local and short (≤4 entries each).
- Hartford trajectory is rendered as a polyline with markers; off-corridor
  points sit on the x-axis at E_PEP=0 with X markers, not at the curve.
- The cliff crossing (Δt_pharm=10h) is annotated explicitly rather than
  with a single from-to arrow.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Constants matching pharmacy_sensitivity.py
T_CRIT_UPPER = 34.0  # ρ=1.0
T_CRIT_LOWER = 32.0  # ρ=0.30

# Cities to label on Panel A.  Kept short to avoid clutter.
HIGHLIGHT_CITIES = {
    "Milwaukee":   dict(dx=+0.6, dy=+0.04, ha="left"),
    "Jackson":     dict(dx=-0.6, dy=+0.06, ha="right"),
    "SanJuan":     dict(dx=+0.8, dy=-0.08, ha="left"),
    "Houston":     dict(dx=+0.8, dy=-0.06, ha="left"),
    "PalmBeachCO": dict(dx=+0.8, dy=+0.08, ha="left"),
    "NewHaven":    dict(dx=-0.8, dy=-0.06, ha="right"),
    "Hartford":    dict(dx=-0.8, dy=+0.10, ha="right"),
}

# Color palette
C_UPPER     = "#1f4e79"   # solid dark blue, upper envelope
C_LOWER     = "#5b8bbf"   # mid blue, lower envelope
C_CORRIDOR  = "#a8c8e4"   # light blue fill
C_CLIFF_HI  = "#a83232"   # dark red, ρ=1.0 cliff
C_CLIFF_LO  = "#cc7a7a"   # rose, ρ=0.30 cliff
C_CITIES    = "#222222"   # near-black scatter
C_HARTFORD  = "#8e2424"   # burgundy trajectory
C_HARTFORD_OFF = "#000000"  # off-corridor markers


def _draw_corridor(ax, corridor_df, *, label_corridor=True):
    """Draw the shaded envelope + bounding curves + t_crit cliffs."""
    ax.fill_between(
        corridor_df["t_acq_h"],
        corridor_df["E_PEP_lower"],
        corridor_df["E_PEP_upper"],
        alpha=0.30, color=C_CORRIDOR,
        label="Envelope corridor" if label_corridor else None,
    )
    ax.plot(corridor_df["t_acq_h"], corridor_df["E_PEP_upper"],
            color=C_UPPER, linewidth=2.0,
            label=r"Upper: $\rho=1.0$")
    ax.plot(corridor_df["t_acq_h"], corridor_df["E_PEP_lower"],
            color=C_LOWER, linewidth=1.4, linestyle="--",
            label=r"Lower: $\rho=0.30$")
    ax.axvline(T_CRIT_UPPER, color=C_CLIFF_HI, linestyle=":", linewidth=1.5,
               alpha=0.85,
               label=fr"$t_{{\rm crit}}(\rho=1.0)={T_CRIT_UPPER:.1f}$h")
    ax.axvline(T_CRIT_LOWER, color=C_CLIFF_LO, linestyle=":", linewidth=1.0,
               alpha=0.6,
               label=fr"$t_{{\rm crit}}(\rho=0.30)={T_CRIT_LOWER:.1f}$h")


def _panel_a(ax, corridor_df, city_positions_df):
    """Panel A: Full corridor + 34 city baselines."""
    _draw_corridor(ax, corridor_df, label_corridor=True)

    baseline = city_positions_df[city_positions_df["delta_t_pharm_h"] == 0]
    on_corridor = baseline[~baseline["past_tcrit"]]

    ax.scatter(on_corridor["t_acq_h"], on_corridor["E_PEP_upper"],
               color=C_CITIES, s=22, alpha=0.75, zorder=3,
               label=f"{len(on_corridor)} US metros at $\\Delta t_{{\\rm pharm}}=0$")

    # Annotate highlight cities with thin tether lines
    for city_name, opts in HIGHLIGHT_CITIES.items():
        row = baseline[baseline["city"] == city_name]
        if row.empty:
            continue
        x = float(row["t_acq_h"].iloc[0])
        y = float(row["E_PEP_upper"].iloc[0])
        ax.annotate(
            city_name,
            xy=(x, y),
            xytext=(x + opts["dx"], y + opts["dy"]),
            fontsize=8.5, color="black",
            ha=opts["ha"], va="center",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.55),
        )

    ax.set_xlabel(
        r"Acquisition time  $t_{\rm acq}$  (hours)",
        fontsize=10.5,
    )
    ax.set_ylabel(r"$E_{\rm PEP}$", fontsize=10.5)
    ax.set_title("A.  Envelope corridor with 34 high-burden US metros at baseline",
                 fontsize=10.5, loc="left", pad=8)
    ax.set_xlim(-0.5, 38)
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.25)

    leg = ax.legend(loc="lower left", fontsize=8.0, framealpha=0.92,
                    handlelength=2.0, labelspacing=0.45)
    leg.get_frame().set_edgecolor("lightgray")


def _panel_b(ax, corridor_df, city_positions_df):
    """Panel B: Hartford pharmacy displacement trajectory (zoomed)."""
    # Background corridor (no labels — Panel A has those)
    _draw_corridor(ax, corridor_df, label_corridor=False)

    hartford = (
        city_positions_df[city_positions_df["city"].str.contains("Hartford", case=False, na=False)]
        .sort_values("delta_t_pharm_h")
        .reset_index(drop=True)
    )

    on_corridor = hartford[~hartford["past_tcrit"]]
    off_corridor = hartford[hartford["past_tcrit"]]

    # Trajectory line connecting on-corridor positions
    if len(on_corridor) >= 2:
        ax.plot(on_corridor["t_acq_h"], on_corridor["E_PEP_upper"],
                color=C_HARTFORD, linewidth=1.7, linestyle="-",
                marker="o", markersize=5.5, zorder=4,
                label="Hartford trajectory (on corridor)")

    # Annotate each Δt_pharm value along Hartford's trajectory
    for _, row in on_corridor.iterrows():
        dtp = int(row["delta_t_pharm_h"])
        x = float(row["t_acq_h"])
        y = float(row["E_PEP_upper"])
        ax.annotate(
            f"{dtp}h",
            xy=(x, y),
            xytext=(x + 0.4, y + 0.025),
            fontsize=7.5, color=C_HARTFORD,
        )

    # Off-corridor markers — sit on the x-axis at y=0
    if not off_corridor.empty:
        ax.scatter(off_corridor["t_acq_h"],
                   np.zeros(len(off_corridor)),
                   color=C_HARTFORD_OFF, marker="X", s=90,
                   edgecolor="white", linewidth=1.0, zorder=5,
                   label=f"Hartford off corridor (Δt_pharm ≥10h)")
        for _, row in off_corridor.iterrows():
            dtp = int(row["delta_t_pharm_h"])
            x = float(row["t_acq_h"])
            ax.annotate(
                f"{dtp}h",
                xy=(x, 0.0),
                xytext=(x + 0.3, -0.025),
                fontsize=7.5, color=C_HARTFORD_OFF,
            )

    # Cliff-crossing call-out: Δt_pharm=10h, the first off-corridor sample
    crossing = hartford[hartford["delta_t_pharm_h"] == 10]
    if not crossing.empty:
        cx = float(crossing["t_acq_h"].iloc[0])
        ax.annotate(
            "Cliff crossing\nat $\\Delta t_{\\rm pharm}=10$h",
            xy=(cx, 0.015),
            xytext=(cx - 6.5, 0.32),
            fontsize=9, color=C_HARTFORD, weight="bold", ha="left",
            arrowprops=dict(arrowstyle="->", color=C_HARTFORD,
                            lw=1.2, alpha=0.85),
        )

    ax.set_xlabel(
        r"Acquisition time  $t_{\rm acq}$  (hours)",
        fontsize=10.5,
    )
    ax.set_ylabel(r"$E_{\rm PEP}$", fontsize=10.5)
    ax.set_title("B.  Hartford pharmacy displacement (Δt_struct=24.4h)",
                 fontsize=10.5, loc="left", pad=8)
    # Zoom into the Hartford range
    ax.set_xlim(22, 38)
    ax.set_ylim(-0.06, 0.62)
    ax.grid(True, alpha=0.25)

    leg = ax.legend(loc="upper right", fontsize=8.0, framealpha=0.92,
                    handlelength=2.0, labelspacing=0.45)
    leg.get_frame().set_edgecolor("lightgray")


def plot_corridor(corridor_df: pd.DataFrame,
                  city_positions_df: pd.DataFrame,
                  out_path: Path) -> None:
    fig, (axL, axR) = plt.subplots(
        nrows=1, ncols=2, figsize=(14, 6),
        gridspec_kw=dict(width_ratios=[1.55, 1.0], wspace=0.22),
    )
    _panel_a(axL, corridor_df, city_positions_df)
    _panel_b(axR, corridor_df, city_positions_df)

    fig.suptitle(
        "Pharmacy delay displaces cities along a fixed kinetic boundary",
        fontsize=12.5, y=1.02,
    )
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(str(out_path).replace(".png", ".pdf"), bbox_inches="tight")
    print(f"Saved: {out_path}")
    print(f"Saved: {str(out_path).replace('.png', '.pdf')}")


def main():
    in_dir = Path("v3_revision/results/pharmacy_sensitivity_corrected")
    corridor_df = pd.read_csv(in_dir / "envelope_corridor.csv")
    city_positions_df = pd.read_csv(in_dir / "city_envelope_positions.csv")

    out_path = in_dir / "Figure_S14_envelope_corridor.png"
    plot_corridor(corridor_df, city_positions_df, out_path)


if __name__ == "__main__":
    main()
