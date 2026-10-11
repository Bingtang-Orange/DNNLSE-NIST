# DNNLSE-NIST

Reference MATLAB implementation of the numerical inverse scattering transform (NIST) for the defocusing nonlocal nonlinear Schrödinger equation (DNNLSE).

## Quick Setup

Before running any script, add all subdirectories to the MATLAB path by running `setup.m` from the repository root:

```matlab
setup
```

---

## Directory Tree

```text
DNNLSE-NIST/
├── .gitignore
├── LICENSE
├── README.md
├── setup.m
├── example/
│   ├── decay_ray.m
│   ├── quickstart_direct_scattering.m
│   ├── rhp1_evolution.m
│   ├── rhp2_evolution.m
│   ├── rhp2_vs_rhp4.m
│   ├── rhp4_spectral_convergence.m
│   ├── rhp4_t30_simulation.m
│   ├── run_fig_reflection_coefficients.m
│   ├── run_fig_step_error_benchmark.m
│   ├── run_fig_truncation_test.m
│   ├── run_step_validation_table.m
│   ├── t30_core_zoom.m
│   └── t30_full_domain.m
├── src/
│   ├── cauchy/
│   │   ├── algorithm_5_4.m
│   │   ├── assemble_cauchy_block.m
│   │   ├── assemble_global.m
│   │   ├── cauchy_def_5_1.m
│   │   ├── cauchy_def_5_3.m
│   │   ├── cauchy_def_5_5.m
│   │   ├── cauchy_def_5_6.m
│   │   ├── cauchy_def_5_7.m
│   │   ├── cauchy_def_5_8.m
│   │   ├── cauchy_def_5_9.m
│   │   ├── cauchy_def_5_10.m
│   │   ├── cauchy_def_5_11.m
│   │   ├── cauchy_def_5_12.m
│   │   ├── cauchy_def_5_13.m
│   │   ├── cauchy_def_5_14.m
│   │   ├── cauchy_kernel_interval.m
│   │   ├── cheby.m
│   │   └── psi_row.m
│   ├── direct_scattering/
│   │   ├── compute_S.m
│   │   └── compute_reflection.m
│   ├── geometry/
│   │   ├── create_curve_def.m
│   │   ├── curve_derivative.m
│   │   ├── curve_mapping.m
│   │   └── forward_map.m
│   └── rhp_solvers/
│       ├── clencurt.m
│       ├── compute_delta.m
│       ├── eval_dnnlse_asymptotics.m
│       ├── solve_rhp1.m
│       ├── solve_rhp2.m
│       └── solve_rhp4.m
└── test_and_data/
    ├── rhp1_vs_rhp4.m
    └── t30_rhp4_data.mat
```
