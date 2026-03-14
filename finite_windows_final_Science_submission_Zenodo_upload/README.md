# Prevention Theorem: Biological and Structural Barriers to Post-Exposure Prophylaxis (PEP)

This repository contains the source code, datasets, and reproduction pipeline for the manuscript:
**"Biological and Structural Barriers to Post-Exposure Prophylaxis: The Prevention Theorem for Parenteral HIV Exposure"** (Science, 2026).

## Abstract
The "Prevention Theorem" establishes a rigorous mathematical framework for the timing and efficacy of Post-Exposure Prophylaxis (PEP). While mucosal HIV exposure (e.g., sexual transmission) provides a multi-day "window of opportunity" for PEP, parenteral exposure (e.g., injection drug use) significantly compresses this window. This repository provides the tools to reproduce our analysis of how source viral load uncertainty and city-level structural delays interact to determine PEP efficacy across the United States.

## Repository Structure
The repository is organized as a "flat" structure for maximum compatibility and ease of use:

- `PEP_parenteral_perelson.py`: Models parenteral PEP efficacy and generates **Figure 1**.
- `PEP_stochastic_perelson.py`: Monte Carlo analysis of source viral load uncertainty and generates **Figure 2**.
- `aidsvu_city_profiles.py`: Processes 34-city AIDSVu datasets to derive structural delays and efficacy.
- `city_stratified_figures.py`: Generates city-level comparisons (**Figure 3** and **Supplemental Figure S3**).
- `Makefile`: Orchestrates the full reproduction pipeline.
- `requirements.txt`: Python dependencies.
- `aidsvu datasets/`: Directory containing city-level epidemiological data.
- `archive/`: Contains legacy scripts, preprints, and intermediate experimental data.

## Installation

### Prerequisites
- Python 3.9 or higher
- `pip` (Python package manager)
- `make` (standard Unix build tool)

### Setup
1. Clone or download the repository.
2. Install the required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

## Reproducing Results

We provide a `Makefile` to simplify the generation of all manuscript figures and datasets.

### Generate All Figures
To regenerate the entire analysis (Figures 1, 2, 3, and Supplemental Figures), run:
```bash
make all
```

### Individual Components
- **Figure 1 (Parenteral VL Sweep):** `make fig1`
- **Figure 2 (Stochastic Uncertainty):** `make fig2`
- **Figure 3 & Supp S3 (City Analysis):** `make fig3`

### Outputs
After running the pipeline, the following files will be generated:
- `pep_parenteral_vl_sweep.png` (Fig 1)
- `pep_stochastic_vl_uncertainty.png` (Fig 2)
- `results/city_analysis/Fig_CityStratified_PEP.png` (Fig 3)
- `results/city_analysis/Fig_CityComparison_Focus.png` (Supp Fig S3)
- `pep_stochastic_results/`: Directory containing detailed regime analysis and numerical data.

## Data Sources
- **AIDSVu**: Metropolitan area data (2023-2025) for viral suppression and HIV prevalence among PWID.
- **CDC/SSP**: Syringe Service Program (SSP) coverage data used for structural delay derivations.

## License
This project is licensed under the terms of the MIT License.

## Citation
If you use this code or data in your research, please cite the manuscript:
*Demidont, A. C. (2026). Biological and Structural Barriers to Post-Exposure Prophylaxis. Science.*
