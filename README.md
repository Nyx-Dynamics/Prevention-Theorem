# The Prevention Theorem

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Time-Dependent Constraints on Post-Exposure Prophylaxis for HIV**

## Overview

This repository contains the computational implementation and analysis code for the Prevention Theorem, which formalizes HIV prevention as a mathematical boundary condition problem. The theorem establishes that post-exposure prophylaxis (PEP) can achieve true prevention (R₀(e) = 0) only when initiated within a finite biological window prior to irreversible proviral integration.

### Key Findings

- **Prevention Condition**: True prevention requires R₀(e) = 0, corresponding to zero probability of productive infection
- **Time-Dependent Efficacy**: PEP efficacy decays monotonically as a function of time post-exposure
- **Mucosal Window**: ~72 hours for sexual exposures (buffered by tissue barriers)
- **Parenteral Window**: ~12-24 hours for injection exposures (bypasses mucosal bottlenecks)
- **Irreversible Transition**: Once proviral integration occurs, the system enters an irreducible infection state

## Repository Structure

```
Prevention-Theorem/
├── SRC/
│   ├── prevention_theorem_figures.py    # Main figure generation
│   ├── PEP_mucosal.py                   # Mucosal vs parenteral analysis
│   ├── generate_prevention_data.py      # Data generation utilities
│   └── data_output/                     # Generated data files
├── Figures/
│   ├── Figure_1_Prevention_Theorem_Dynamics.tif
│   └── Figure_2_Window_Compression.tif
├── preprints_org/
│   ├── prevention_theorem_preprints.tex # Manuscript (LaTeX)
│   └── prevention_theorem_preprints.pdf # Compiled manuscript
├── graphical_abstracts.py               # Graphical abstract generation
└── README.md
```

## Mathematical Framework

### The Prevention Theorem

For any viral exposure *e*, true prevention is defined as:

```
R₀(e) = 0
```

This condition implies zero probability of establishing a productive, transmissible infection.

### Time-Dependent Efficacy

PEP efficacy is modeled as a function of cumulative biological transitions:

```
E_PEP(t) = (1 - P_seed(t))·ε_max + (P_seed(t) - P_int(t))·ε_mid + P_int(t)·ε_min
```

Where:
- `P_seed(t)` = Cumulative probability of reservoir seeding
- `P_int(t)` = Cumulative probability of proviral integration
- `ε_max`, `ε_mid`, `ε_min` = Efficacy coefficients by phase

### PEP Window Corollary

Post-exposure prophylaxis can enforce R₀(e) = 0 if and only if initiated within a finite biological window `t < t_crit` prior to irreversible proviral integration.

## Installation

```bash
git clone https://github.com/Nyx-Dynamics/Prevention-Theorem.git
cd Prevention-Theorem
pip install -r requirements.txt
```

### Dependencies

- Python ≥ 3.9
- NumPy ≥ 1.21
- SciPy ≥ 1.7
- Matplotlib ≥ 3.5
- Pandas ≥ 1.3
- Seaborn ≥ 0.11

## Usage

### Generate Figures

```python
cd SRC
python prevention_theorem_figures.py
```

### Generate Graphical Abstracts

```python
python graphical_abstracts.py
```

## Citation

If you use this code or the Prevention Theorem framework in your research, please cite:

### Paper Citation

```bibtex
@article{demidont2026prevention,
  author = {Demidont, A.C.},
  title = {The Prevention Theorem: Time-Dependent Constraints on Post-Exposure Prophylaxis for HIV},
  journal = {Preprints.org},
  year = {2026},
  doi = {10.20944/preprints202601.XXXX.v1},
  url = {https://www.preprints.org/manuscript/202601.XXXX}
}
```

### Software Citation

```bibtex
@software{demidont2026prevention_code,
  author = {Demidont, A.C.},
  title = {Prevention-Theorem: Computational Implementation of Time-Dependent HIV Prevention Constraints},
  year = {2026},
  publisher = {GitHub},
  version = {v1.0.0},
  url = {https://github.com/Nyx-Dynamics/Prevention-Theorem},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

### Software Dependencies

```bibtex
@article{harris2020numpy,
  author = {Harris, Charles R. and others},
  title = {Array programming with NumPy},
  journal = {Nature},
  year = {2020},
  volume = {585},
  pages = {357--362},
  doi = {10.1038/s41586-020-2649-2}
}

@article{virtanen2020scipy,
  author = {Virtanen, Pauli and others},
  title = {SciPy 1.0: fundamental algorithms for scientific computing in Python},
  journal = {Nature Methods},
  year = {2020},
  volume = {17},
  pages = {261--272},
  doi = {10.1038/s41592-019-0686-2}
}

@article{hunter2007matplotlib,
  author = {Hunter, John D.},
  title = {Matplotlib: A 2D graphics environment},
  journal = {Computing in Science \& Engineering},
  year = {2007},
  volume = {9},
  number = {3},
  pages = {90--95},
  doi = {10.1109/MCSE.2007.55}
}
```

## Related Work

- **PWID HIV Prevention Analysis**: [github.com/Nyx-Dynamics/HIV_Prevention_PWID](https://github.com/Nyx-Dynamics/HIV_Prevention_PWID)
- **Algorithmic Bias Epidemiology**: [github.com/Nyx-Dynamics/algorithmic-bias-epidemiology-academic](https://github.com/Nyx-Dynamics/algorithmic-bias-epidemiology-academic)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**A.C. Demidont, DO**
Independent Researcher
Nyx Dynamics, LLC
Email: ac.demidont@nyxdynamics.com

## Acknowledgments

The author thanks the HIV prevention research community whose published work informed model parameterization, and the people who inject drugs (PWID) community advocates whose testimony informed characterization of structural barriers.

---

*This research was conducted independently and released as open-source work. The author reports prior employment with Gilead Sciences, Inc. (2020-2024); Gilead had no role in this research.*
