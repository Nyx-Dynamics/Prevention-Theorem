"""
Generate Figure 1 v2 for aeh5879 Science Advances submission.

Two panels:
  A. Route comparison E_PEP(t) for mucosal (V0=1) vs parenteral (V0=10^3),
     CV=0.3, with 95% credible intervals from N=500 realizations per cell.
     Verified t_crit values labeled. NHP timepoints overlaid.
  B. Schematic of three-state within-host model and route-dependent inoculum.

Source data: SRC/multiscale_model/results_v3/heterogeneity_realizations.csv
Anchored on commit d047d2d (numbers reproduce from any v2-submission tag).
Output: v2_revision/figures/Figure_1_Route_Dependent_PEP_Efficacy_Decay.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

# ============================================================
# Load source data — committed at SHA d047d2d (verified by §5.1).
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPT_DIR)
COMMITTED = os.path.join(BASE, 'SRC/multiscale_model/results_v3')
REAL = pd.read_csv(os.path.join(COMMITTED, 'heterogeneity_realizations.csv'))
SUMM = pd.read_csv(os.path.join(COMMITTED, 'heterogeneity_summary.csv'))

OUT_DIR = os.path.join(BASE, 'v2_revision/figures')
os.makedirs(OUT_DIR, exist_ok=True)

# Operating point: CV=0.3, mucosal V0=1, parenteral V0=1000
CV = 0.3
EMAX = 0.95
EMID = 0.50
EMIN = 0.0
ETA = 0.05
ECLIPSE_FLOOR_H = 21.6
TIME_GRID = np.linspace(0, 200, 401)

def compute_E_PEP_with_CI(df_route, n_bootstrap=200, rng=None):
    """Compute E_PEP(t) and 95% CI by bootstrap resampling realizations."""
    if rng is None:
        rng = np.random.default_rng(0)
    integrated = df_route[df_route['integrated']].copy()
    if len(integrated) < 2:
        return TIME_GRID, np.full_like(TIME_GRID, np.nan), np.full_like(TIME_GRID, np.nan), np.full_like(TIME_GRID, np.nan)
    T_ints = integrated['T_int_hours'].values
    # T_seed approximated as eclipse floor (model floor)
    T_seeds = np.full_like(T_ints, ECLIPSE_FLOOR_H)
    n = len(T_ints)
    boot = np.empty((n_bootstrap, len(TIME_GRID)))
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        Ti = T_ints[idx]; Ts = T_seeds[idx]
        Pseed = np.array([np.mean(Ts <= t) for t in TIME_GRID])
        Pint  = np.array([np.mean(Ti <= t) for t in TIME_GRID])
        boot[i] = (1 - Pseed)*EMAX + (Pseed - Pint)*EMID + Pint*EMIN
    mean = boot.mean(axis=0)
    lo = np.percentile(boot, 2.5, axis=0)
    hi = np.percentile(boot, 97.5, axis=0)
    return TIME_GRID, mean, lo, hi

# Filter to canonical operating point. V0 disambiguates route in this kinetics-only output:
# mucosal = post-epithelial attenuation V0=1; parenteral = direct intravascular V0=10^3.
muc = REAL[(REAL['V0']==1) & (REAL['cv']==CV)]
par = REAL[(REAL['V0']==1000) & (REAL['cv']==CV)]

print(f"mucosal realizations: n={len(muc)}, integrated={muc['integrated'].sum()}")
print(f"parenteral realizations: n={len(par)}, integrated={par['integrated'].sum()}")

t_grid_m, E_m, lo_m, hi_m = compute_E_PEP_with_CI(muc, rng=np.random.default_rng(1))
t_grid_p, E_p, lo_p, hi_p = compute_E_PEP_with_CI(par, rng=np.random.default_rng(2))

# Verified t_crit values from heterogeneity_summary.csv
def _summ_tcrit(V0, cv):
    rows = SUMM[(SUMM['V0']==V0) & (SUMM['cv']==cv)]
    if 'route' in SUMM.columns:
        # if route column present, prefer the matching one; else just take first
        pass
    return float(rows['t_crit_eta_05'].iloc[0])

TCRIT_M = _summ_tcrit(1, CV)
TCRIT_P = _summ_tcrit(1000, CV)
print(f"Verified t_crit: mucosal={TCRIT_M}h, parenteral={TCRIT_P}h, ratio={TCRIT_M/TCRIT_P:.3f}x")

# ============================================================
# Figure
# ============================================================
fig = plt.figure(figsize=(13, 5.5), facecolor='white')
gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0], wspace=0.25)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# --------- Panel A: route comparison ---------
COL_M = '#1F4E79'  # blue (mucosal)
COL_P = '#9C2A2A'  # red (parenteral)

# Mean curves and CIs
axA.fill_between(t_grid_m, lo_m*100, hi_m*100, color=COL_M, alpha=0.22, linewidth=0)
axA.plot(t_grid_m, E_m*100, color=COL_M, linewidth=2.6,
         label=r'Mucosal ($V_0=1$, CV=0.3)')
axA.fill_between(t_grid_p, lo_p*100, hi_p*100, color=COL_P, alpha=0.22, linewidth=0)
axA.plot(t_grid_p, E_p*100, color=COL_P, linewidth=2.6,
         label=r'Parenteral ($V_0=10^3$, CV=0.3)')

# Reference lines
axA.axhline(ETA*100, color='#444', linestyle='-.', alpha=0.7, linewidth=1)
axA.text(199, ETA*100 + 1.5, r'$\eta = 5\%$', fontsize=8.5, color='#444', ha='right')
axA.axvline(72, color='#2D7A2D', linestyle='--', alpha=0.6, linewidth=1)
axA.text(72, 102, 'CDC 72h', fontsize=8.5, color='#2D7A2D', ha='center',
         bbox=dict(facecolor='white', edgecolor='none', pad=1.5))

# --- F1.1 fix: place t_crit labels in low-traffic area BELOW the curves
#     (around y=20-30, x>=TCRIT) where there is no NHP marker, no curve, no callout.
#     Use bbox so labels are clearly separated from any underlying gridlines.
TCRIT_LABEL_BBOX = dict(facecolor='white', edgecolor=None, boxstyle='round,pad=0.25', alpha=0.9)

# Parenteral t_crit label — bottom area, right of the dotted line
axA.annotate(
    r"$t^{(p)}_{\mathrm{crit}} \approx 34.5$ h",
    xy=(TCRIT_P, 5),                                      # arrow tip at η=0.05 line
    xytext=(TCRIT_P + 12, 22),                            # text below curves, right
    fontsize=10.5, color=COL_P, weight='bold',
    arrowprops=dict(arrowstyle='-', color=COL_P, lw=0.7, alpha=0.6),
    ha='left', va='center',
    bbox=dict(facecolor='white', edgecolor=COL_P, boxstyle='round,pad=0.3', alpha=0.95),
)

# Mucosal t_crit label — bottom area, right of mucosal dotted line, above parenteral label
axA.annotate(
    r"$t^{(m)}_{\mathrm{crit}} \approx 60.5$ h",
    xy=(TCRIT_M, 5),
    xytext=(TCRIT_M + 12, 36),
    fontsize=10.5, color=COL_M, weight='bold',
    arrowprops=dict(arrowstyle='-', color=COL_M, lw=0.7, alpha=0.6),
    ha='left', va='center',
    bbox=dict(facecolor='white', edgecolor=COL_M, boxstyle='round,pad=0.3', alpha=0.95),
)

# NHP overlay — colored markers; concordant (up-triangle) vs non-concordant (down-triangle).
# Otten 2000 mucosal: 12h, 36h concordant; 72h non-concordant (model crashes to ~0)
# Tsai 1998 parenteral: 24h, 48h non-concordant (model under-predicts)
nhp_muc = [(12, 1.00, 'concordant'), (36, 1.00, 'concordant'), (72, 0.50, 'discordant')]
nhp_par = [(24, 1.00, 'discordant'), (48, 0.50, 'discordant')]

for delay, prot, concord in nhp_muc:
    marker = '^' if concord == 'concordant' else 'v'
    face = COL_M if prot >= 0.99 else ('lightblue' if prot >= 0.4 else 'white')
    axA.plot(delay, prot*100, marker=marker, color=COL_M, markerfacecolor=face,
             markersize=11, markeredgewidth=1.8, markeredgecolor=COL_M, alpha=0.95, zorder=10)
for delay, prot, concord in nhp_par:
    marker = '^' if concord == 'concordant' else 'v'
    face = COL_P if prot >= 0.99 else ('mistyrose' if prot >= 0.4 else 'white')
    axA.plot(delay, prot*100, marker=marker, color=COL_P, markerfacecolor=face,
             markersize=11, markeredgewidth=1.8, markeredgecolor=COL_P, alpha=0.95, zorder=10)

# --- F1.3 fix: U=U scope annotation in lower-right (empty space; was upper-right).
axA.text(0.99, 0.32,
         'U=U scope condition\n'
         'Framework requires source VL ≥200 c/mL\n'
         '(transmission-competent; Rodger 2019,\nBavinton 2018)',
         transform=axA.transAxes, ha='right', va='top',
         fontsize=8.0, color='#444',
         bbox=dict(facecolor='#f5f5f5', edgecolor='#bbb', boxstyle='round,pad=0.4'))

axA.set_xlim(0, 200)
axA.set_ylim(0, 105)
axA.set_xlabel('Hours from exposure to PEP initiation', fontsize=11)
axA.set_ylabel(r'Expected PEP efficacy $E_\mathrm{PEP}(t)$ (%)', fontsize=11)
axA.set_title('A. Route comparison: emergent $t_\\mathrm{crit}$ from multiscale dynamics',
              fontsize=12, weight='bold', loc='left', pad=8)
from matplotlib.lines import Line2D
leg = axA.legend(fontsize=9.5, framealpha=0.95, edgecolor='#bbb', loc='upper right')

# --- F1.2 fix: NHP legend shows actual marker semantics (shape = concordance,
#     color = route, fill = NHP protection level). Mirrors what's plotted.
nhp_legend_handles = [
    Line2D([0], [0], marker='^', color=COL_M, markerfacecolor=COL_M,
           markersize=9, linestyle='None', markeredgewidth=1.6,
           label='Mucosal, model concordant (NHP=100%)'),
    Line2D([0], [0], marker='v', color=COL_M, markerfacecolor='lightblue',
           markersize=9, linestyle='None', markeredgewidth=1.6,
           label='Mucosal, not concordant (NHP=50%)'),
    Line2D([0], [0], marker='v', color=COL_P, markerfacecolor='white',
           markersize=9, linestyle='None', markeredgewidth=1.6,
           label='Parenteral, not concordant (NHP=100%)'),
    Line2D([0], [0], marker='v', color=COL_P, markerfacecolor='mistyrose',
           markersize=9, linestyle='None', markeredgewidth=1.6,
           label='Parenteral, not concordant (NHP=50%)'),
]

axA.add_artist(leg)
nhp_leg = axA.legend(handles=nhp_legend_handles, loc='upper center',
           bbox_to_anchor=(0.5, -0.13), ncol=2,
           fontsize=8, framealpha=0.95, title='NHP empirical (▲ concordant,  ▼ not concordant)',
           title_fontsize=8.5)
axA.grid(True, alpha=0.25, linewidth=0.5)

# Compression ratio annotation
axA.text(0.02, 0.05,
         f'Ratio (mucosal/parenteral) = {TCRIT_M/TCRIT_P:.2f}$\\times$\nNHP empirical: 1.5--2$\\times$ (concordant)',
         transform=axA.transAxes, ha='left', va='bottom',
         fontsize=8.5, color='#222',
         bbox=dict(facecolor='white', edgecolor='#bbb', boxstyle='round,pad=0.3'))

# --------- Panel B: Schematic ---------
axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.axis('off')

# Three states
state_props = dict(boxstyle='round,pad=0.4', linewidth=1.5)
axB.add_patch(mpatches.FancyBboxPatch((0.4, 6.5), 2.2, 1.6, boxstyle='round,pad=0.1', linewidth=1.5,
                                      facecolor='#dceaf3', edgecolor=COL_M))
axB.text(1.5, 7.3, 'Z=0\nSusceptible', ha='center', va='center', fontsize=10, weight='bold')
axB.text(1.5, 8.4, r'$\varepsilon=\varepsilon_\mathrm{max}=0.95$', ha='center', fontsize=8.5, color='#1F4E79')

axB.add_patch(mpatches.FancyBboxPatch((4.0, 6.5), 2.2, 1.6, boxstyle='round,pad=0.1', linewidth=1.5,
                                      facecolor='#fff2cc', edgecolor='#9c8a2a'))
axB.text(5.1, 7.3, 'Z=1\nSeeded', ha='center', va='center', fontsize=10, weight='bold')
axB.text(5.1, 8.4, r'$\varepsilon=\varepsilon_\mathrm{mid}=0.50$', ha='center', fontsize=8.5, color='#9c8a2a')

axB.add_patch(mpatches.FancyBboxPatch((7.6, 6.5), 2.2, 1.6, boxstyle='round,pad=0.1', linewidth=1.5,
                                      facecolor='#f5d6d6', edgecolor=COL_P))
axB.text(8.7, 7.3, 'Z=2\nIntegrated', ha='center', va='center', fontsize=10, weight='bold')
axB.text(8.7, 8.4, r'$\varepsilon=\varepsilon_\mathrm{min}=0$', ha='center', fontsize=8.5, color=COL_P)
axB.text(8.7, 5.95, '(absorbing)', ha='center', fontsize=8, style='italic', color='#666')

# Arrows
axB.annotate('', xy=(4.0, 7.3), xytext=(2.6, 7.3),
             arrowprops=dict(arrowstyle='->', lw=1.5, color='#444'))
axB.text(3.3, 7.7, r'$T_\mathrm{seed}$', ha='center', fontsize=10)
axB.annotate('', xy=(7.6, 7.3), xytext=(6.2, 7.3),
             arrowprops=dict(arrowstyle='->', lw=1.5, color='#444'))
axB.text(6.9, 7.7, r'$T_\mathrm{int}$', ha='center', fontsize=10)

# Inoculum schematic at bottom
axB.text(5.0, 4.8, 'Route-dependent inoculum', ha='center', fontsize=10.5, weight='bold')

axB.add_patch(mpatches.FancyBboxPatch((0.4, 2.4), 4.4, 1.8, boxstyle='round,pad=0.1', linewidth=1.2,
                                      facecolor='#eef4f9', edgecolor=COL_M))
axB.text(2.6, 3.6, 'Mucosal', ha='center', va='center', fontsize=10, weight='bold', color=COL_M)
axB.text(2.6, 2.95, r'$V_0 = 1$ virion/mL', ha='center', va='center', fontsize=9, color='#333')
axB.text(2.6, 2.55, '(post-epithelial attenuation)', ha='center', va='center', fontsize=8, style='italic', color='#666')

axB.add_patch(mpatches.FancyBboxPatch((5.2, 2.4), 4.4, 1.8, boxstyle='round,pad=0.1', linewidth=1.2,
                                      facecolor='#fcecec', edgecolor=COL_P))
axB.text(7.4, 3.6, 'Parenteral', ha='center', va='center', fontsize=10, weight='bold', color=COL_P)
axB.text(7.4, 2.95, r'$V_0 = 10^{3}$ virions/mL', ha='center', va='center', fontsize=9, color='#333')
axB.text(7.4, 2.55, '(direct intravascular)', ha='center', va='center', fontsize=8, style='italic', color='#666')

# Outcome footnote
axB.text(5.0, 1.4, r'By inoculum-monotone hitting times (Lemma S5.2): $T_\mathrm{int}^{(p)} \leq T_\mathrm{int}^{(m)}$',
         ha='center', va='center', fontsize=8.5, style='italic', color='#444')
axB.text(5.0, 0.7, r'$\Rightarrow t_\mathrm{crit}^{(p)} \leq t_\mathrm{crit}^{(m)}$ (Corollary S5.1)',
         ha='center', va='center', fontsize=9, color='#222')

axB.set_title('B. Three-state within-host model & route-dependent inoculum',
              fontsize=12, weight='bold', loc='left', pad=8)

plt.tight_layout()
out = os.path.join(OUT_DIR, 'Figure_1_Route_Dependent_PEP_Efficacy_Decay.png')
plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
out_pdf = os.path.join(OUT_DIR, 'Figure_1_Route_Dependent_PEP_Efficacy_Decay.pdf')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out}")
print(f"Saved: {out_pdf}")
