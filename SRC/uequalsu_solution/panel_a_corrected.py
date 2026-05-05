"""
Corrected Panel A: parenteral PEP efficacy with U=U-respecting framing.

Generates a 2x2 figure replacing the original Panel A:

  A. Conditional PEP efficacy by source VL — parenteral (existing model,
     unchanged). This is what the paper currently shows.

  B. Joint absolute prevention by source VL — parenteral. This is
     p_transmit(VL, parenteral) * P(PEP succeeds | transmit). At VL<50,
     this collapses smoothly toward zero — but never reaches zero,
     because parenteral U=U is not clinically established (Baggaley 2006
     is a single-point-estimate, no PARTNER analogue exists).

  C. Joint absolute prevention by source VL — mucosal (receptive vaginal).
     PARTNER2 floor at VL<200 forces this to exactly zero, recovering
     the U=U statement: no transmission to prevent at suppressed VL.

  D. Side-by-side asymmetry headline at VL=50: sexual joint prevention
     is structurally zero; parenteral is small but non-zero. This is
     the figure that defends the paper's parenteral claims while
     showing the sexual U=U story is respected.

Reads no command-line args. Writes:
    SRC/uequalsu_solution/panel_a_corrected.png
    SRC/uequalsu_solution/panel_a_corrected.pdf
"""

import os
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'route_models'))
sys.path.insert(0, HERE)

from parenteral_route import ParenteralExposureModel
from mucosal_route import InfectionEstablishmentModel
from transmission_probability import (
    p_transmit, joint_prevention, SEXUAL_VL_FLOOR,
)


VL_GRID = [50, 200, 1000, 10_000, 50_000, 100_000, 500_000, 1_000_000]
HOURS_GRID = np.linspace(0, 120, 240)


def parenteral_curves():
    out = {}
    for vl in VL_GRID:
        m = ParenteralExposureModel(source_viral_load=vl,
                                    exposure_type='pwid_shared_needle')
        eff = np.array([m.pep_efficacy(h)['pep_efficacy'] for h in HOURS_GRID])
        out[vl] = eff
    return out


def mucosal_curves():
    out = {}
    m = InfectionEstablishmentModel()
    eff = np.array([m.pep_efficacy(h)['pep_efficacy'] for h in HOURS_GRID])
    for vl in VL_GRID:
        out[vl] = eff
    return out


def vl_label(vl):
    if vl < 100:
        return f'VL {int(vl)} (suppressed)'
    if vl < 1_000:
        return f'VL {int(vl)}'
    if vl < 1e6:
        return f'VL {int(vl/1000)}K'
    return f'VL {vl/1e6:.1f}M (acute)'


def plot():
    par = parenteral_curves()
    muc = mucosal_curves()

    cmap = plt.cm.RdYlGn_r
    norm = mcolors.LogNorm(vmin=50, vmax=1e6)

    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle(
        'Panel A — Corrected: PEP Efficacy with U=U-Respecting Framing\n'
        'Parenteral (no U=U floor) vs Mucosal (PARTNER2 floor at VL<200)',
        fontsize=13, fontweight='bold', y=0.995,
    )

    # ----- A: parenteral conditional efficacy -----
    ax = axes[0, 0]
    for vl in VL_GRID:
        ax.plot(HOURS_GRID, par[vl] * 100,
                color=cmap(norm(vl)), linewidth=2.2, label=vl_label(vl))
    ax.axvline(72, color='k', ls='--', lw=1.2, alpha=0.6)
    ax.axhline(50, color='gray', ls='-.', lw=1, alpha=0.6)
    ax.set_xlabel('Hours from exposure to PEP initiation')
    ax.set_ylabel('Conditional PEP efficacy (%)')
    ax.set_title('A. Parenteral — conditional efficacy\n'
                 'P(PEP succeeds | transmission occurred)',
                 loc='left', fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, 120); ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)

    # ----- B: parenteral joint absolute prevention -----
    ax = axes[0, 1]
    for vl in VL_GRID:
        joint = joint_prevention(par[vl], vl, 'pwid_shared_needle')
        ax.plot(HOURS_GRID, joint * 1000,
                color=cmap(norm(vl)), linewidth=2.2, label=vl_label(vl))
    ax.set_xlabel('Hours from exposure to PEP initiation')
    ax.set_ylabel('Infections prevented per 1,000 exposures')
    ax.set_title('B. Parenteral — joint absolute prevention\n'
                 'p_transmit(VL) x P(PEP succeeds) [no U=U floor]',
                 loc='left', fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, 120)
    ax.grid(True, alpha=0.3)

    # ----- C: mucosal joint absolute prevention with U=U floor -----
    ax = axes[1, 0]
    for vl in VL_GRID:
        joint = joint_prevention(muc[vl], vl, 'receptive_vaginal')
        ax.plot(HOURS_GRID, joint * 1000,
                color=cmap(norm(vl)), linewidth=2.2, label=vl_label(vl))
    ax.text(60, ax.get_ylim()[1] * 0.85,
            f'VL < {SEXUAL_VL_FLOOR}: joint = 0\n'
            '(PARTNER2: zero linked transmissions)',
            fontsize=10, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow',
                      edgecolor='goldenrod'))
    ax.set_xlabel('Hours from exposure to PEP initiation')
    ax.set_ylabel('Infections prevented per 1,000 exposures')
    ax.set_title('C. Mucosal (receptive vaginal) — joint prevention\n'
                 'p_transmit(VL) x P(PEP succeeds) [U=U floor at VL<200]',
                 loc='left', fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, 120)
    ax.grid(True, alpha=0.3)

    # ----- D: asymmetry headline at VL=50 -----
    ax = axes[1, 1]
    routes_demo = [
        ('Parenteral (PWID needle share)', par[50],
         'pwid_shared_needle', 'crimson'),
        ('Mucosal (receptive vaginal)', muc[50],
         'receptive_vaginal', 'steelblue'),
        ('Mucosal (receptive anal)', muc[50],
         'receptive_anal', 'darkblue'),
    ]
    for label, eff_curve, route, color in routes_demo:
        joint = joint_prevention(eff_curve, 50, route)
        ax.plot(HOURS_GRID, joint * 1000,
                color=color, linewidth=2.5, label=label)
    ax.set_xlabel('Hours from exposure to PEP initiation')
    ax.set_ylabel('Infections prevented per 1,000 exposures')
    ax.set_title('D. Asymmetry at VL=50 (suppressed)\n'
                 'Sexual = 0 (U=U).  Parenteral != 0 (gap in evidence).',
                 loc='left', fontweight='bold', fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim(0, 120)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    out_png = os.path.join(HERE, 'panel_a_corrected.png')
    out_pdf = os.path.join(HERE, 'panel_a_corrected.pdf')
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"Saved: {os.path.relpath(out_png, PROJECT_ROOT)}")
    print(f"Saved: {os.path.relpath(out_pdf, PROJECT_ROOT)}")


def print_summary_table():
    print()
    print("Joint prevention at hours_to_pep = 24h, by VL and route")
    print("(infections prevented per 1,000 exposures)")
    print("=" * 72)
    par = parenteral_curves()
    muc = mucosal_curves()
    idx_24h = int(np.argmin(np.abs(HOURS_GRID - 24)))

    print(f"{'VL':>10}  {'Parenteral':>14}  {'Mucosal (RV)':>14}  "
          f"{'Mucosal (RA)':>14}")
    print("-" * 60)
    for vl in VL_GRID:
        j_par = joint_prevention(par[vl][idx_24h], vl,
                                  'pwid_shared_needle') * 1000
        j_rv = joint_prevention(muc[vl][idx_24h], vl,
                                 'receptive_vaginal') * 1000
        j_ra = joint_prevention(muc[vl][idx_24h], vl,
                                 'receptive_anal') * 1000
        print(f"{vl:>10.0f}  {j_par:>14.3f}  {j_rv:>14.3f}  {j_ra:>14.3f}")
    print()
    print("Note row VL=50: mucosal entries are exactly 0 (U=U floor);")
    print("parenteral entry is small but nonzero (gap in evidence).")


if __name__ == '__main__':
    print("=" * 72)
    print("PANEL A CORRECTED — U=U-respecting framing")
    print("=" * 72)
    plot()
    print_summary_table()
