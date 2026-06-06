"""
generate_S2_tables.py — supplement S2 LaTeX table generator

Reads the three CSVs produced by pharmacy_sensitivity.py (corrected v4) and
writes three LaTeX table files for inclusion in the supplement:

  table_S2a.tex  — Envelope corridor at key t_acq values
  table_S2b.tex  — 34 high-burden US metros, baseline positions on corridor
  table_S2c.tex  — Pharmacy displacement summary (cities past t_crit per Δt_pharm)

Each output is a self-contained \\begin{table}...\\end{table} block.
Include in the supplement with \\input{table_S2a} etc.

Run after pharmacy_sensitivity.py (or as part of reproduce_all_v4.py).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
IN_DIR = REPO_ROOT / 'v3_revision' / 'results' / 'pharmacy_sensitivity_corrected'
OUT_DIR = REPO_ROOT / 'v3_revision' / 'tables'
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Table S2a: Envelope corridor at key t_acq values
# ---------------------------------------------------------------------------
def make_table_S2a(corridor_df: pd.DataFrame) -> str:
    """Sample the corridor at clinically/mechanistically interpretable timepoints."""
    sample_t = [0, 6, 12, 18, 24, 28, 30, 32, 33, 34, 35, 36]
    rows = []
    for t in sample_t:
        idx = (corridor_df['t_acq_h'] - t).abs().idxmin()
        r = corridor_df.iloc[idx]
        rows.append({
            't': float(r['t_acq_h']),
            'upper': float(r['E_PEP_upper']),
            'lower': float(r['E_PEP_lower']),
            'width': float(r['envelope_width']),
            'past_upper': bool(r['past_tcrit_upper']),
            'past_lower': bool(r['past_tcrit_lower']),
        })

    body = []
    for r in rows:
        flag = ''
        if r['past_upper'] and r['past_lower']:
            flag = r'\dag\dag'  # past both cliffs
        elif r['past_upper']:
            flag = r'\dag'
        elif r['past_lower']:
            flag = r'\ddag'
        body.append(
            f"  {r['t']:5.1f} & {r['upper']:.4f} & {r['lower']:.4f} & "
            f"{r['width']:.4f} & {flag} \\\\"
        )

    return r"""\begin{table}[!htb]
\centering
\caption{\textbf{Table S2a. Envelope corridor at key acquisition times.}
The corridor is bounded above by $E_{\rm PEP}$ at perfect adherence
($\rho=1.0$, cliff at $t_{\rm crit}=34.0$\,h) and below by low adherence
($\rho=0.30$, cliff at $t_{\rm crit}=32.0$\,h). Reproduced from
\texttt{envelope\_corridor.csv}, sampled at twelve timepoints.
\dag~past $t_{\rm crit}(\rho=1.0)$.
\ddag~past $t_{\rm crit}(\rho=0.30)$.}
\label{tab:S2a_corridor}
\begin{tabular}{rrrrl}
\toprule
$t_{\rm acq}$ (h) & $E_{\rm PEP}^{\rm upper}$ & $E_{\rm PEP}^{\rm lower}$ & Width & Cliff \\
\midrule
""" + '\n'.join(body) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Table S2b: City baseline positions on corridor
# ---------------------------------------------------------------------------
def make_table_S2b(city_positions_df: pd.DataFrame) -> str:
    """All 34 cities at Δt_pharm = 0, sorted by structural delay (ascending)."""
    baseline = (
        city_positions_df[city_positions_df['delta_t_pharm_h'] == 0]
        .sort_values('delta_t_struct_h')
        .reset_index(drop=True)
    )

    # Split into two side-by-side columns for compactness (17 rows each)
    half = (len(baseline) + 1) // 2
    left = baseline.iloc[:half].reset_index(drop=True)
    right = baseline.iloc[half:].reset_index(drop=True)

    rows = []
    for i in range(half):
        L = left.iloc[i]
        L_city = str(L['city']).replace('_', r'\_')
        L_state = ' '.join(str(L.get('state', '')).split())
        L_str = f"{L_city}, {L_state}" if L_state else L_city
        L_text = (
            f"{L_str} & {L['delta_t_struct_h']:.1f} & "
            f"{L['E_PEP_upper']:.3f} & {L['E_PEP_lower']:.3f}"
        )
        if i < len(right):
            R = right.iloc[i]
            R_city = str(R['city']).replace('_', r'\_')
            R_state = ' '.join(str(R.get('state', '')).split())
            R_str = f"{R_city}, {R_state}" if R_state else R_city
            R_text = (
                f"{R_str} & {R['delta_t_struct_h']:.1f} & "
                f"{R['E_PEP_upper']:.3f} & {R['E_PEP_lower']:.3f}"
            )
        else:
            R_text = " &  &  & "
        rows.append(f"  {L_text} & {R_text} \\\\")

    return r"""\begin{table}[!htb]
\centering
\small
\caption{\textbf{Table S2b. Baseline corridor positions for 34 high-burden US
metropolitan areas.} Each city is shown at $\Delta t_{\rm pharm} = 0$, so
$t_{\rm acq} = \Delta t_{\rm struct}$. Cities sorted by structural delay
(ascending). $E_{\rm PEP}^{\rm upper}$ is at $\rho=1.0$ and
$E_{\rm PEP}^{\rm lower}$ is at $\rho=0.30$. Reproduced from
\texttt{city\_envelope\_positions.csv}. All 34 cities are on-corridor at
$\Delta t_{\rm pharm}=0$; Hartford is the only city that crosses
$t_{\rm crit}(\rho=1.0)=34.0$\,h within the $\Delta t_{\rm pharm} \in [0,12]$\,h
sweep (see Table~\ref{tab:S2c_displacement}).}
\label{tab:S2b_baseline}
\begin{tabular}{lrrr@{\hspace{1.2em}}lrrr}
\toprule
City & $\Delta t_{\rm struct}$ & $E^{\rm upper}$ & $E^{\rm lower}$ &
City & $\Delta t_{\rm struct}$ & $E^{\rm upper}$ & $E^{\rm lower}$ \\
 & (h) & & & & (h) & & \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Table S2c: Pharmacy displacement summary
# ---------------------------------------------------------------------------
def make_table_S2c(summary_df: pd.DataFrame) -> str:
    """The pharmacy delay sweep: cities displaced past t_crit at each Δt_pharm."""
    rows = []
    for _, r in summary_df.iterrows():
        displaced = r['displaced_city_names']
        if displaced == '-' or pd.isna(displaced):
            displaced = '---'
        rows.append(
            f"  {int(r['delta_t_pharm_h']):>2}h & "
            f"{int(r['n_cities_displaced_past_tcrit'])}/34 & "
            f"{r['median_t_acq_h']:.1f} & "
            f"{r['mean_E_PEP_upper']:.3f} & "
            f"{r['mean_E_PEP_lower']:.3f} & "
            f"{displaced} \\\\"
        )

    return r"""\begin{table}[!htb]
\centering
\caption{\textbf{Table S2c. Pharmacy displacement summary.} For each
hypothetical pharmacy dispensing delay $\Delta t_{\rm pharm}$, the table
reports the number of cities displaced past $t_{\rm crit}(\rho=1.0) = 34.0$\,h
(off corridor), the cohort median $t_{\rm acq}$, and the cohort means of
upper and lower envelope positions. Cohort means are unweighted averages
of city-level kinetic predictions; cities displaced past $t_{\rm crit}$
contribute zero. These are not population-weighted preventability estimates.
Reproduced from \texttt{pharmacy\_displacement\_summary.csv}.}
\label{tab:S2c_displacement}
\begin{tabular}{rrrrrl}
\toprule
$\Delta t_{\rm pharm}$ & Off corridor & Median $t_{\rm acq}$ &
Cohort mean & Cohort mean & Displaced \\
 & ($n$ of 34) & (h) & $E^{\rm upper}$ & $E^{\rm lower}$ & city/cities \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    corridor_df = pd.read_csv(IN_DIR / 'envelope_corridor.csv')
    city_positions_df = pd.read_csv(IN_DIR / 'city_envelope_positions.csv')
    summary_df = pd.read_csv(IN_DIR / 'pharmacy_displacement_summary.csv')

    artifacts = [
        ('table_S2a.tex', make_table_S2a(corridor_df)),
        ('table_S2b.tex', make_table_S2b(city_positions_df)),
        ('table_S2c.tex', make_table_S2c(summary_df)),
    ]
    for fname, content in artifacts:
        out_path = OUT_DIR / fname
        out_path.write_text(content)
        print(f"Wrote {out_path}")

    print(f"\nInclude in supplement with:")
    print(f"  \\input{{v3_revision/tables/table_S2a}}")
    print(f"  \\input{{v3_revision/tables/table_S2b}}")
    print(f"  \\input{{v3_revision/tables/table_S2c}}")


if __name__ == '__main__':
    main()
