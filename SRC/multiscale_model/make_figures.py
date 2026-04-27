"""
Generate figures and validation outputs for the multiscale model.

Figures produced:
  Fig A: Efficacy decay E_PEP(t) by V0 (route compression, derived not assumed)
  Fig B: T_int distributions (showing extinction at low V0, tight distributions
         at high V0)
  Fig C: t_crit values vs log10(V0) — the compression curve
  Fig D: Compare derived t_crit to manuscript-claimed and NHP data
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

mpl.rcParams['font.family'] = 'DejaVu Sans'

# Load results
summary = pd.read_csv('/home/claude/multiscale_ode/results/mc_summary.csv')
realizations = pd.read_csv('/home/claude/multiscale_ode/results/mc_realizations.csv')
curves = pd.read_csv('/home/claude/multiscale_ode/results/mc_efficacy_curves.csv')

V0_grid = sorted(summary['V0'].unique())

# Color gradient by V0 (blue = mucosal, red = parenteral acute)
cmap = plt.cm.RdYlBu_r
colors = {V0: cmap(np.log10(V0)/4.0) for V0 in V0_grid}

fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
fig.suptitle("Multiscale Within-Host Model: t_crit Emerges from Dynamics",
             fontsize=14, weight='bold', y=0.995)

# ===== Panel A: E_PEP(t) curves =====
ax = axes[0, 0]
for V0 in V0_grid:
    col = f'E_PEP_V0_{V0}'
    if col in curves.columns:
        label = (f'V₀ = 1 (mucosal)' if V0 == 1
                 else f'V₀ = 10³ (parenteral)' if V0 == 1000
                 else f'V₀ = 10⁴ (acute)' if V0 == 10000
                 else f'V₀ = {V0}')
        ax.plot(curves['time_h'], curves[col]*100,
                color=colors[V0], linewidth=2.0,
                label=label if V0 in [1, 10, 100, 1000, 10000] else None)
ax.axhline(5, color='#666', linestyle=':', alpha=0.7)
ax.text(195, 6, 'η = 5%', fontsize=8, color='#666', ha='right')
ax.axvline(72, color='#444', linestyle=':', alpha=0.4)
ax.text(72, 100, 'CDC 72h', fontsize=8, color='#444', ha='center',
        bbox=dict(facecolor='white', edgecolor='none', pad=1))
ax.set_xlim(0, 200)
ax.set_ylim(0, 105)
ax.set_xlabel('Hours from exposure to PEP initiation', fontsize=10.5)
ax.set_ylabel('Expected PEP efficacy (%)', fontsize=10.5)
ax.set_title('A. Efficacy decay across V₀\n(curves emerge from ODE+stochastic dynamics)',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=9, framealpha=0.9, edgecolor='#ccc')
ax.grid(True, alpha=0.25, linewidth=0.5)

# ===== Panel B: T_int distributions =====
ax = axes[0, 1]
for V0 in [1, 10, 100, 1000, 10000]:
    sub = realizations[(realizations['V0']==V0) & realizations['integrated']]
    if len(sub) > 5:
        ax.hist(sub['T_int_hours'], bins=30, alpha=0.55,
                color=colors[V0],
                label=f'V₀ = {V0:>5} (n={len(sub)})')
ax.set_xlabel('T_int (hours)', fontsize=10.5)
ax.set_ylabel('Count', fontsize=10.5)
ax.set_title('B. Distribution of integration completion times\n'
             '(T_int = first time R(t) ≥ 10 cells/mL)',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.25, linewidth=0.5)

# ===== Panel C: t_crit vs log10 V0 =====
ax = axes[1, 0]
log_V0 = summary['log10_V0'].values

for eta, color, lw in [(0.80, '#3D7A3D', 2),
                        (0.50, '#1F4E79', 2.5),
                        (0.10, '#B8860B', 2),
                        (0.05, '#9C2A2A', 2.5)]:
    col = f't_crit_eta_{int(eta*100):02d}'
    ax.plot(log_V0, summary[col], 'o-', color=color, linewidth=lw,
            markersize=7, markerfacecolor='white', markeredgewidth=2,
            label=f'η = {eta}')

# Annotations
ax.axhline(72, color='#888', linestyle=':', alpha=0.5)
ax.text(0, 73, 'CDC 72h', fontsize=8.5, color='#666', ha='left')

# Manuscript claim band
ax.axhspan(16, 28, xmin=0.7, xmax=0.85, alpha=0.15, color='red')
ax.text(3.5, 22, 'Manuscript\nparenteral\nclaim: 16–28h',
        fontsize=8, color='#9C2A2A', ha='center', va='center')
ax.axhspan(68, 76, xmin=0, xmax=0.15, alpha=0.15, color='blue')
ax.text(0, 72, 'Manuscript\nmucosal\nclaim: 68–76h',
        fontsize=8, color='#1F4E79', ha='center', va='center')

ax.set_xlabel('log₁₀ V₀ (virions/mL at exposure)', fontsize=10.5)
ax.set_ylabel('t_crit (hours)', fontsize=10.5)
ax.set_title('C. Critical-time-to-PEP-failure as a function of route\n'
             'Derived from multiscale dynamics, not hand-set parameters',
             fontsize=11, weight='bold', loc='left')
ax.legend(loc='upper right', fontsize=9.5, framealpha=0.95)
ax.grid(True, alpha=0.25, linewidth=0.5)
ax.set_xlim(-0.3, 4.3)
ax.set_ylim(0, 90)

# ===== Panel D: extinction probability + handoff time =====
ax = axes[1, 1]
ax2 = ax.twinx()

ax.plot(log_V0, summary['p_extinct']*100, 'o-', color='#9C2A2A',
        linewidth=2.5, markersize=8, markerfacecolor='white',
        markeredgewidth=2, label='P(extinction)')
ax.set_xlabel('log₁₀ V₀ (virions/mL at exposure)', fontsize=10.5)
ax.set_ylabel('Extinction probability (%)', fontsize=10.5, color='#9C2A2A')
ax.tick_params(axis='y', labelcolor='#9C2A2A')
ax.set_xlim(-0.3, 4.3)
ax.set_ylim(0, 100)

ax2.plot(log_V0, summary['T_seed_median_h'], 's-', color='#1F4E79',
         linewidth=2.5, markersize=8, markerfacecolor='white',
         markeredgewidth=2, label='T_seed median (handoff time)')
ax2.set_ylabel('T_seed median (hours)', fontsize=10.5, color='#1F4E79')
ax2.tick_params(axis='y', labelcolor='#1F4E79')
ax2.axhline(21.6, color='#888', linestyle=':', alpha=0.5)
ax2.text(4.0, 22.3, 'eclipse phase floor (Perelson 1996)',
         fontsize=7.5, color='#666', ha='right', style='italic')
ax2.set_ylim(0, 60)

ax.set_title('D. The two mechanisms of route compression\n'
             '(founder-bottleneck delay vs. eclipse-phase floor)',
             fontsize=11, weight='bold', loc='left')
ax.grid(True, alpha=0.25, linewidth=0.5)

# Combined legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labels1+labels2, loc='center right', fontsize=9)

plt.tight_layout()
plt.savefig('/home/claude/multiscale_ode/results/multiscale_results.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('/home/claude/multiscale_ode/results/multiscale_results.pdf',
            bbox_inches='tight', facecolor='white')
plt.close()

print("Saved: results/multiscale_results.png/.pdf")

# Print a markdown-formatted comparison table
print("\n" + "="*72)
print("COMPARISON: published phenomenological model vs. new multiscale ODE")
print("="*72)

# Manuscript claims vs new model
comparisons = [
    ('Mucosal (V₀=1)',         '68-76h',         summary[summary['V0']==1]['t_crit_eta_05'].iloc[0]),
    ('Parenteral PWID (V₀=10³)', '16-28h',         summary[summary['V0']==1000]['t_crit_eta_05'].iloc[0]),
    ('Parenteral acute (V₀=10⁴)', '~12-16h',         summary[summary['V0']==10000]['t_crit_eta_05'].iloc[0]),
]
print(f"{'Route':<28} {'Manuscript claim':<20} {'New ODE-derived':<20}")
print("-"*68)
for r, m, n in comparisons:
    print(f"{r:<28} {m:<20} {n:.1f}h")

m_ratio = summary[summary['V0']==1]['t_crit_eta_05'].iloc[0] / \
          summary[summary['V0']==1000]['t_crit_eta_05'].iloc[0]
print(f"\nCompression ratio (mucosal/parenteral): {m_ratio:.2f}x")
print(f"Manuscript-claimed compression: 3.0x")
