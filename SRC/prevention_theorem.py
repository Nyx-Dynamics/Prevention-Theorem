import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --- EPIDEMICS JOURNAL STYLING ---
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8


def save_fig(fig, name):
    fig.savefig(f"{name}.eps", format='eps', bbox_inches='tight')
    fig.savefig(f"{name}.tif", format='tif', dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
    print(f"Saved {name}")


# --- FIGURE 1: THEORETICAL FRAMEWORK ---
def plot_figure_1():
    t = np.linspace(0, 96, 500)  # Hours

    # Sigmoidal functions for biological probabilities
    # Seed: Reservoir seeding starts early
    p_seed = 1 / (1 + np.exp(-0.15 * (t - 24)))
    # Int: Integration lags seeding
    p_int = 1 / (1 + np.exp(-0.2 * (t - 48)))

    # Efficacy Function E_PEP(t)
    # E_max = 0.95, E_mid = 0.5, E_min = 0
    e_max, e_mid, e_min = 0.95, 0.5, 0.0
    e_pep = (1 - p_seed) * e_max + (p_seed - p_int) * e_mid + p_int * e_min

    # Irreversibility Metric (1 - E_PEP)
    irreversibility = 1 - e_pep

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])

    # Panel A: Establishment Timeline
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t, p_seed, 'b--', label='Reservoir Seeding $P_{seed}(t)$')
    ax1.plot(t, p_int, 'r-', linewidth=2, label='Proviral Integration $P_{int}(t)$')
    ax1.set_ylabel('Cumulative Probability')
    ax1.set_title('A. Infection Establishment Dynamics')
    ax1.legend(loc='lower right', frameon=False, fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Panel B: Time-Dependent Efficacy
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(t, e_pep, 'g-', linewidth=2.5)
    ax2.axvline(72, color='gray', linestyle=':')
    ax2.text(74, 0.8, 'Critical Window\n(~72h Mucosal)', fontsize=9, color='gray')
    ax2.set_ylabel('PEP Efficacy $E_{PEP}(t)$')
    ax2.set_title('B. Time-Dependent Prevention Operator')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    # Panel C: Phase Transition to Irreducibility
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.fill_between(t, 0, irreversibility, color='firebrick', alpha=0.2)
    ax3.plot(t, irreversibility, 'k-', linewidth=1.5)
    ax3.set_xlabel('Time Post-Exposure (Hours)')
    ax3.set_ylabel('Probability of Irreversible Infection')
    ax3.set_title('C. Transition to Irreducible State')
    ax3.text(10, 0.1, 'Preventable State\n$R_0(e)=0$ Feasible', color='green')
    ax3.text(60, 0.8, 'Irreducible State\n$R_0(e)>0$ Fixed', color='firebrick')

    # Panel D: Theorem Text Visualization
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    text_str = (
            r"$\bf{The\ Prevention\ Theorem}$" + "\n\n" +
            r"Prevention requires $R_0(e,t) = 0$." + "\n\n" +
            r"$R_0(e,t) = 1 - E_{PEP}(t)$" + "\n\n" +
            r"As $t \to \infty$, $P_{int}(t) \to 1 \Rightarrow E_{PEP}(t) \to 0$." + "\n\n" +
            r"$\therefore$ Prevention is only possible" + "\n" +
            r"while $P_{int}(t) < 1$."
    )
    ax4.text(0.5, 0.5, text_str, ha='center', va='center', fontsize=12,
             bbox=dict(boxstyle="round,pad=1", fc="white", ec="black", alpha=0.8))
    ax4.set_title('D. Mathematical Formalism')

    plt.tight_layout()
    save_fig(fig, 'Figure_1_Prevention_Theorem_Dynamics')
    plt.close()


# --- FIGURE 2: MUCOSAL VS PARENTERAL ---
def plot_figure_2():
    t = np.linspace(0, 96, 500)

    # Efficacy decay
    # Mucosal: Slow decay (immune bottlenecks)
    eff_mucosal = 1 / (1 + np.exp(0.15 * (t - 60)))
    # Parenteral: Fast decay (systemic access)
    eff_parenteral = 1 / (1 + np.exp(0.3 * (t - 18)))

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t, eff_mucosal, 'b-', linewidth=2, label='Mucosal Exposure (Sexual)')
    ax.plot(t, eff_parenteral, 'r--', linewidth=2, label='Parenteral Exposure (Injection)')

    # Thresholds
    ax.axvline(24, color='red', linestyle=':', alpha=0.5)
    ax.axvline(72, color='blue', linestyle=':', alpha=0.5)

    ax.text(25, 0.6, 'Parenteral\nWindow\n~24h', color='red', fontsize=9)
    ax.text(73, 0.6, 'Mucosal\nWindow\n~72h', color='blue', fontsize=9)

    ax.set_xlabel('Time from Exposure to PEP (Hours)')
    ax.set_ylabel('Probability of Prevention Success')
    ax.set_title('Figure 2: Compression of Prevention Window for Parenteral Exposure')
    ax.set_xlim(0, 96)
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    save_fig(fig, 'Figure_2_Window_Compression')
    plt.close()


# --- GRAPHICAL ABSTRACT: BOUNDARY-CONDITION TIMELINE ---
def plot_graphical_abstract():
    fig, ax = plt.subplots(figsize=(10, 4))

    # Time axis
    t_max = 100
    ax.axhline(0.5, color='black', linewidth=2)
    ax.set_xlim(-5, t_max + 5)
    ax.set_ylim(0, 1)

    # Key time points
    points = {
        0: ("Exposure", "green"),
        24: ("Parenteral\nBoundary", "red"),
        72: ("Mucosal\nBoundary", "blue"),
        96: ("Irreducible\nState", "black")
    }

    for t, (label, color) in points.items():
        ax.plot(t, 0.5, 'o', color=color, markersize=10)
        ax.text(t, 0.55, label, ha='center', va='bottom', color=color, fontweight='bold')

    # Arrows/Zones
    ax.annotate('', xy=(24, 0.45), xytext=(0, 0.45),
                arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
    ax.text(12, 0.4, 'High Risk\nWindow', ha='center', color='red', fontsize=9)

    ax.annotate('', xy=(72, 0.35), xytext=(0, 0.35),
                arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
    ax.text(36, 0.3, 'Mucosal Window of Opportunity', ha='center', color='blue', fontsize=9)

    # Formalism
    ax.text(50, 0.8, r"$\mathcal{P} \text{ Theorem: } E_{PEP}(t) \to 0 \text{ as } t \to \tau_{critical}$",
            ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    ax.axis('off')
    save_fig(fig, 'Graphical_Abstract_Boundary_Timeline')
    plt.close()


if __name__ == "__main__":
    plot_figure_1()
    plot_figure_2()
    plot_graphical_abstract()
