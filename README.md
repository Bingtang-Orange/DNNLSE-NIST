# DNNLSE-NIST

Reference MATLAB implementation of the high-order numerical inverse scattering transform (NIST) for the defocusing nonlocal nonlinear Schrödinger equation (DNNLSE).

## Quick Setup

Before running any script, initialize all paths from the repository root:
```matlab
setup


DNNLSE-NIST/
├── setup.m                     # Automated MATLAB path initializer
├── LICENSE                     # MIT License
├── README.md                   # Repository file index
├── .gitignore                  # Git ignore rules
├── src/                        # Core numerical solver libraries
│   ├── cauchy/                 # Discrete Cauchy singular operators (Olver 2012)
│   ├── direct_scattering/      # Lax pair ODE solver via Chebfun
│   ├── geometry/               # Conformal contour mappings and differential metrics
│   └── rhp_solvers/            # Matrix RHP solvers and asymptotic evaluators
├── example/                    # Numerical experiments and paper reproduction scripts
└── test_and_data/              # Unit verification tests and precomputed datasets
