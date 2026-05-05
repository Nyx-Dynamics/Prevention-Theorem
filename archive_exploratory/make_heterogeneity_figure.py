"""
Figure: heterogeneity widens T_int distributions, route compression ratio
1.5-1.8x emerges from dynamics (matches NHP, not manuscript's 3x claim).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'

import sys
sys.path.insert(0, '/home/claude/multiscale_v3')
from multiscale_v3 import WithinHostParameters
from run_targeted_heterogeneity import run_mc_targeted, E_PEP_at, t_crit_at, assess


def E_PEP_curve(df, time_grid):
    return np.array([E_PEP_at(df, t) for t in time_grid])


# Two heterogeneity profiles:
#   Profile P0: noise-free (CV=0 everywhere) — matches v2 baseline
#   Profile P1: best-effort heterogeneity (Scenario C from sweep)
#               alpha CV=0.7, T0 CV=0.5, tau_eclipse CV=0.3, others CV=0.3

profiles = {
    'P0_noisefree': {},
    'P1_full_het': {
        'beta': 0.3, 'c': 0.3, 'delta': 0.3,
        'alpha': 0.7, 'T0': 0.5, 'tau_eclipse': 0.3,
    },
}

V0_grid = [1, 1000]
N = 800

# Run all combinations
results = {}
for prof_name, cv_dict in profiles.items():
    for V0 in V0_grid:
        print(f"Running {prof_name}, V0={V0}...")
        df = run_mc_targeted(V0, cv_dict, N)
        results[(prof_name, V0)] = df

# Build figure
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
fig.suptitle("Multiscale Model with Inter-Individual Heterogeneity",
             fontsize=14, weight='bold', y=0.995)

time_grid = np.linspace(0, 200, 401)

# Panel A: T_int distributions
ax = axes[0, 0]
for prof_name, color, label in [('P0_noisefree', '#888888', 'No heterogeneity'),
                                  ('P1_full_het', '#1F4E79', 'Full heterogeneity')]:
    for V0, alpha_val, ls in [(1, 0.45, '-'), (1000, 0.65, '-')]:
        df = results[(prof_name, V0)]
        integrated = df[df['integrated']]
        if len(integrated) > 5:
            label_full = f"V₀={V0}, {label}" if V0 == 1000 else None
            ax.hist(integrated['T_int_hours'], bins=40, alpha=alpha_val,
                    color=color if V0 == 1000 else 'lightcoral' if 'P0' in prof_name else '#9C2A2A',
                    label=f"V₀={V0}, {label}",
                    histtype='stepfilled', edgecolor='black', linewidth=0.4)
ax.set_xlabel('T_int (hours)', fontsize=10.5)
ax.set_ylabel('Count (N=800 per cell)', fontsize=10.5)
ax.set_title('A. T_int distributions widen with heterogeneity\n'
             '(but median timing is preserved)',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=9, framealpha=0.92)
ax.grid(True, alpha=0.25, linewidth=0.5)
ax.set_xlim(0, 120)

# Panel B: E_PEP curves
ax = axes[0, 1]
for V0, base_color in [(1, '#1F4E79'), (1000, '#9C2A2A')]:
    for prof_name, ls, alpha_val in [('P0_noisefree', '--', 0.8),
                                       ('P1_full_het', '-', 1.0)]:
        df = results[(prof_name, V0)]
        eff = E_PEP_curve(df, time_grid)
        label = (f"V₀={V0}, "
                 f"{'noise-free' if 'P0' in prof_name else 'heterogeneous'}")
        ax.plot(time_grid, eff*100, ls, color=base_color,
                linewidth=2.0, alpha=alpha_val, label=label)

ax.axhline(5, color='#666', linestyle=':', alpha=0.6)
ax.text(195, 6.5, 'η=5%', fontsize=8, color='#666', ha='right')
ax.axvline(72, color='#444', linestyle=':', alpha=0.4)
ax.text(72, 100, 'CDC 72h', fontsize=8, color='#444', ha='center',
        bbox=dict(facecolor='white', edgecolor='none', pad=1))
# NHP timepoints
nhp_pts = [(24, 1.00, 'Tsai IV', '#9C2A2A'),
           (48, 0.50, 'Tsai IV', '#9C2A2A'),
           (12, 1.00, 'Otten', '#1F4E79'),
           (36, 1.00, 'Otten', '#1F4E79'),
           (72, 0.50, 'Otten', '#1F4E79')]
for delay, target, study, color in nhp_pts:
    ax.scatter(delay, target*100, s=110, marker='v',
               edgecolor=color, facecolor='white', linewidth=1.6, zorder=5)
ax.scatter([], [], s=110, marker='v', edgecolor='black',
           facecolor='white', linewidth=1.6, label='NHP protection data')

ax.set_xlim(0, 200); ax.set_ylim(0, 105)
ax.set_xlabel('Hours from exposure to PEP initiation', fontsize=10.5)
ax.set_ylabel('Expected PEP efficacy (%)', fontsize=10.5)
ax.set_title('B. Efficacy curves: route + heterogeneity\n'
             'NHP triangles overlaid (Tsai 1998, Otten 2000)',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=8.5, framealpha=0.92, loc='upper right')
ax.grid(True, alpha=0.25, linewidth=0.5)

# Panel C: Compression ratio across scenarios
ax = axes[1, 0]
scenarios = []
for prof_name in ['P0_noisefree', 'P1_full_het']:
    df_m = results[(prof_name, 1)]
    df_p = results[(prof_name, 1000)]
    t_05_m = t_crit_at(df_m, 0.05)
    t_05_p = t_crit_at(df_p, 0.05)
    t_50_m = t_crit_at(df_m, 0.50)
    t_50_p = t_crit_at(df_p, 0.50)
    scenarios.append({
        'profile': prof_name,
        'mucosal_t_crit_05': t_05_m,
        'parenteral_t_crit_05': t_05_p,
        'compression_05': t_05_m / t_05_p if t_05_p > 0 else float('nan'),
        'compression_50': t_50_m / t_50_p if t_50_p > 0 else float('nan'),
    })

# Bar plot of compression ratios
labels = ['Noise-free\n(v2 baseline)', 'Heterogeneous\n(v3, full)']
comp_05 = [s['compression_05'] for s in scenarios]
comp_50 = [s['compression_50'] for s in scenarios]
xpos = np.arange(len(labels))
w = 0.35
b1 = ax.bar(xpos - w/2, comp_05, w, color='#9C2A2A', alpha=0.85,
             label='at η=0.05', edgecolor='black', linewidth=0.6)
b2 = ax.bar(xpos + w/2, comp_50, w, color='#1F4E79', alpha=0.85,
             label='at η=0.50', edgecolor='black', linewidth=0.6)

# Reference lines
ax.axhline(3.0, color='#9C2A2A', linestyle=':', alpha=0.6, linewidth=1.5)
ax.text(1.5, 3.05, 'Manuscript claim: 3.0×',
        fontsize=8.5, color='#9C2A2A', ha='right', style='italic')
ax.axhline(1.5, color='#3D7A3D', linestyle='--', alpha=0.6, linewidth=1.5)
ax.text(1.5, 1.55, 'NHP-observed: ~1.5× (Tsai/Otten)',
        fontsize=8.5, color='#3D7A3D', ha='right', style='italic')

# Bar labels
for bars, vals in [(b1, comp_05), (b2, comp_50)]:
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x()+bar.get_width()/2, v+0.05,
                f'{v:.2f}×', ha='center', fontsize=9, weight='bold')

ax.set_xticks(xpos); ax.set_xticklabels(labels)
ax.set_ylabel('Mucosal/parenteral t_crit ratio', fontsize=10.5)
ax.set_ylim(0, 3.5)
ax.set_title('C. Route compression ratio: derived vs. claimed\n'
             'Multiscale model gives ~1.6-1.8×, matching NHP, not 3×',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.25, linewidth=0.5, axis='y')

# Panel D: NHP concordance comparison
ax = axes[1, 1]
nhp_labels_plot = []
nhp_targets = []
nhp_p0_model = []
nhp_p1_model = []

nhp_eval_pts = [
    ('Tsai IV, 24h', 'parenteral', 24, 1.00),
    ('Tsai IV, 48h', 'parenteral', 48, 0.50),
    ('Otten, 12h',   'mucosal',    12, 1.00),
    ('Otten, 36h',   'mucosal',    36, 1.00),
    ('Otten, 72h',   'mucosal',    72, 0.50),
]
for name, route, delay, target in nhp_eval_pts:
    V0 = 1 if route == 'mucosal' else 1000
    df_p0 = results[('P0_noisefree', V0)]
    df_p1 = results[('P1_full_het', V0)]
    nhp_labels_plot.append(name)
    nhp_targets.append(target)
    nhp_p0_model.append(E_PEP_at(df_p0, delay))
    nhp_p1_model.append(E_PEP_at(df_p1, delay))

xpos = np.arange(len(nhp_labels_plot))
w = 0.27
ax.bar(xpos - w, nhp_targets, w, color='#3D7A3D', alpha=0.85,
       label='NHP observed', edgecolor='black', linewidth=0.6)
ax.bar(xpos, nhp_p0_model, w, color='#aaaaaa', alpha=0.85,
       label='Model (noise-free)', edgecolor='black', linewidth=0.6)
ax.bar(xpos + w, nhp_p1_model, w, color='#1F4E79', alpha=0.85,
       label='Model (heterogeneous)', edgecolor='black', linewidth=0.6)

ax.set_xticks(xpos)
ax.set_xticklabels(nhp_labels_plot, rotation=20, ha='right', fontsize=9)
ax.set_ylabel('Protection / E_PEP', fontsize=10.5)
ax.set_ylim(0, 1.15)
ax.set_title('D. NHP concordance by timepoint\n'
             'Heterogeneity narrows but does not close the gap',
             fontsize=11, weight='bold', loc='left')
ax.legend(fontsize=9.5, framealpha=0.92)
ax.grid(True, alpha=0.25, linewidth=0.5, axis='y')

plt.tight_layout()
plt.savefig('/home/claude/multiscale_v3/heterogeneity_comparison.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('/home/claude/multiscale_v3/heterogeneity_comparison.pdf',
            bbox_inches='tight', facecolor='white')
plt.close()

print("\nSaved heterogeneity_comparison.png/.pdf")
print("\nCompression ratio summary:")
for s in scenarios:
    print(f"  {s['profile']:<18}: t_crit ratio at η=0.05 = {s['compression_05']:.2f}x, "
          f"at η=0.50 = {s['compression_50']:.2f}x")

# Save a summary CSV for later use
sumrows = []
for prof_name in ['P0_noisefree', 'P1_full_het']:
    for V0 in V0_grid:
        df = results[(prof_name, V0)]
        a = assess(df, 'mucosal' if V0 == 1 else 'parenteral')
        sumrows.append({
            'profile': prof_name, 'V0': V0,
            'p_extinct': a['p_extinct'], 'n_int': a['n_int'],
            'T_int_median': a['T_int_median'],
            'T_int_p5': a['T_int_p5'], 'T_int_p95': a['T_int_p95'],
            't_crit_05': a['t_crit_05'], 't_crit_50': a['t_crit_50'],
        })
sumdf = pd.DataFrame(sumrows)
sumdf.to_csv('/home/claude/multiscale_v3/heterogeneity_summary.csv', index=False)
print(f"\nSaved heterogeneity_summary.csv\n")
print(sumdf.to_string(index=False))
