# DNNLSE-NIST: High-Order Numerical Inverse Scattering for the Defocusing Nonlocal NLS Equation

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository provides the reference MATLAB implementation of the high-order numerical inverse scattering transform (NIST) for the defocusing nonlocal nonlinear Schrödinger equation (DNNLSE).

The framework implements an operator-theoretic spectral collocation solver based on Olver's approach for Riemann--Hilbert problems (RHPs). It incorporates dynamic steepest-descent scaling, compactified conformal coordinates, and zero-sum compatibility regularizations.

---

## Directory Structure

```text
DNNLSE-NIST/
├── README.md                              # Repository overview and reproduction guide
├── LICENSE                                # MIT License
├── setup.m                                # Automated MATLAB path initializer
├── src/                                   # Core numerical libraries
│   ├── direct_scattering/                 # Lax pair ODE solver via Chebfun
│   │   ├── compute_reflection.m           # High-level reflection coefficients interface
│   │   └── compute_S.m                    # Collocation solver for 2x2 scattering matrix S
│   ├── geometry/                          # Conformal mappings and contour calculus
│   │   ├── create_curve_def.m             # Geometry specification (segments, rays, arcs)
│   │   ├── curve_derivative.m             # Conformal mapping derivatives & metric factors
│   │   ├── curve_mapping.m                # Inverse conformal maps z = M^{-1}(x)
│   │   └── forward_map.m                  # Forward conformal maps x = M(z)
│   ├── cauchy/                            # Discrete Cauchy singular operators (Olver 2012)
│   │   ├── cauchy_def_5_1.m to 5_14.m     # Pairwise boundary Cauchy operators
│   │   ├── cauchy_kernel_interval.m       # Discrete C^+ on reference interval [-1, 1]
│   │   ├── cheby.m                        # Chebyshev polynomial evaluations
│   │   ├── psi_row.m                      # Exterior Cauchy basis evaluations
│   │   ├── algorithm_5_4.m                # Singular endpoint correction vectors
│   │   ├── assemble_cauchy_block.m        # Pairwise block dispatcher (Olver Table 1)
│   │   └── assemble_global.m              # Global block singular Cauchy assembler
│   └── rhp_solvers/                       # Matrix RHP collocation and potential recovery
│       ├── compute_delta.m                # Scalar RH solver for cut jump delta(lambda)
│       ├── clencurt.m                     # Clenshaw--Curtis quadrature weights
│       ├── solve_rhp1.m                   # RHP I solver (short-time physical recovery)
│       ├── solve_rhp2.m                   # RHP II solver (5-contour steepest descent)
│       └── solve_rhp4.m                   # RHP IV solver (9-contour dynamic lens)
├── tests/                                 # Unit verification tests
│   └── test_scalar_rhp_delta.m            # Accuracy verification for scalar delta(lambda)
└── examples/                              # Complete paper reproduction scripts
    ├── quickstart_direct_scattering.m     # Single-point direct scattering
    ├── run_step_validation_table.m        # Closed-form step benchmark Table
    ├── run_fig_step_error_benchmark.m     # Error distribution across 4 step profiles
    ├── run_fig_truncation_test.m          # Truncation domain length L convergence
    ├── run_fig_reflection_coefficients.m  # Real/imaginary/modulus reflection profiles
    ├── run_fig_rhp1_evolution.m           # Short-time wavepacket profiles (t = 0, 1, 2)
    ├── run_fig_rhp2_evolution.m           # Multi-scale dynamics and envelopes
    ├── run_fig_rhp4_t100_envelope.m       # Large-time (t = 100) profile & zero-sum test
    ├── run_fig_rhp4_t500_demod.m          # Demodulated slow-varying profile at t = 500
    ├── run_fig_cross_validation.m         # Cross-stage validation (RHP II vs RHP IV)
    ├── run_fig_p_convergence.m            # Geometric p-spectral convergence (10^-12)
    └── run_fig_decay_slope_benchmark.m    # x = -3t anomalous decay rate fitting
