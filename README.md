# Finite Prevention Windows Under Irreversible Infection Establishment

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Mathematical Framework for HIV Post-Exposure Prophylaxis Timing**

## Overview

This repository contains the computational implementation and analysis code for the Prevention Theorem, which formalizes HIV prevention as a mathematical boundary condition problem. The theorem establishes that post-exposure prophylaxis (PEP) can achieve true prevention (R₀(e) = 0) only when initiated within a finite biological window prior to irreversible proviral integration.

### Key Findings

- **Prevention Condition**: True prevention requires R₀(e) = 0, corresponding to zero probability of productive infection
- **Time-Dependent Efficacy**: PEP efficacy decays monotonically as a function of time post-exposure
- **Mucosal Window**: ~72 hours for sexual exposures (buffered by tissue barriers)
- **Parenteral Window**: ~12–24 hours for injection exposures (bypasses mucosal bottlenecks)
- **Irreversible Transition**: Once proviral integration occurs, the system enters an irreducible infection state

## Repository Structure

```
Prevention-Theorem/
├── SRC/
│   ├── prevention_theorem_figures.py    # Theoretical framework figures (Fig 1–2)
│   ├── route_compression.py            # JID manuscript figures (route-specific windows)
│   ├── PEP_mucosal.py                  # Mucosal vs parenteral PEP analysis
│   ├── middle_ground.py                # Continuous inoculum spectrum model
│   ├── prevention_theorem.py           # Core theorem implementation
│   ├── jid_figures.py                  # JID-specific figure generation
│   ├── generate_prevention_data.py     # Data generation utilities
│   └── graphical_abstracts.py          # Graphical abstract generation
├── figures/                            # All generated figures (PNG, TIFF, EPS)
├── data/                               # Generated data outputs
├── preprints_org/                      # Preprint manuscript (LaTeX/PDF)
├── jid_submission/                     # JID submission materials
├── reproduce_all.py                    # Reproducibility suite (single command)
├── requirements.txt
├── CITATION.cff
├── LICENSE
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

## Reproducibility

All results can be reproduced from a fresh clone:

```bash
git clone https://github.com/Nyx-Dynamics/Prevention-Theorem.git
cd Prevention-Theorem
pip install -r requirements.txt
python reproduce_all.py
```

The reproducibility suite runs all analyses, generates all figures, and validates key manuscript claims with PASS/FAIL checks.

### Generate Individual Components

```bash
# Theoretical framework figures (Fig 1–2)
cd SRC && python prevention_theorem_figures.py

# JID manuscript figures (route-specific windows)
cd SRC && python route_compression.py

# PEP mucosal vs parenteral analysis
cd SRC && python PEP_mucosal.py

# Continuous inoculum spectrum model
cd SRC && python middle_ground.py
```

### Dependencies

- Python ≥ 3.9
- NumPy ≥ 1.21
- SciPy ≥ 1.7
- Matplotlib ≥ 3.5
- Pandas ≥ 1.3
- Seaborn ≥ 0.11

## Citation

### Paper

```bibtex
@article{demidont2026prevention,
  author = {Demidont, A.C.},
  title = {Finite Prevention Windows Under Irreversible Infection Establishment:
           A Mathematical Framework for HIV Post-Exposure Prophylaxis Timing},
  journal = {Preprints.org},
  year = {2026},
  doi = {10.20944/preprints202601.1090.v1},
  url = {https://www.preprints.org/manuscript/202601.1090}
}
```

### Software

```bibtex
@software{demidont2026prevention_code,
  author = {Demidont, A.C.},
  title = {Prevention-Theorem: Computational Implementation of
           Time-Dependent HIV Prevention Constraints},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/Nyx-Dynamics/Prevention-Theorem}
}
```

## Related Work

- **PWID Structural Barriers**: [Nyx-Dynamics/HIV_Prevention_PWID](https://github.com/Nyx-Dynamics/HIV_Prevention_PWID)
- **Algorithmic Bias Epidemiology**: [Nyx-Dynamics/algorithmic-bias-epidemiology-academic](https://github.com/Nyx-Dynamics/algorithmic-bias-epidemiology-academic)
- **LAI-PrEP Bridge Tool**: [Nyx-Dynamics/lai-prep-bridge-tool-pub](https://github.com/Nyx-Dynamics/lai-prep-bridge-tool-pub)

## License

MIT License — see [LICENSE](LICENSE) for details.

## Author

**A.C. Demidont, DO**
Nyx Dynamics, LLC
Email: acdemidont@nyxdynamics.org

---

*This research was conducted independently. The author reports prior employment with Gilead Sciences, Inc. (2020–2024); Gilead had no role in this research.*
