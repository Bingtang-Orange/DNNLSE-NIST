function C_Gamma = cauchy_def_5_6(x, n, def)
%CAUCHY_DEF_5_6 Compute discrete Cauchy operator C^+ on bounded curves.
%   C_GAMMA = CAUCHY_DEF_5_6(X, N, DEF) computes the modified discrete
%   Cauchy operator matrix on a bounded curve.
%
%   Inputs:
%       x       - Collocation nodes on the reference unit interval
%       n       - Number of Chebyshev collocation nodes
%       def     - Bounded curve definition structure
%
%   Outputs:
%       C_Gamma - Discrete Cauchy operator matrix on the curve

    [C_II, ~, ~, ~, ~] = cauchy_kernel_interval(x, n);
    [~, log_dM_L, log_dM_R] = curve_derivative(def, x);

    C_Gamma = C_II;
    C_Gamma(1, 1)     = C_Gamma(1, 1)     - (log_dM_L / (2*1i*pi));
    C_Gamma(end, end) = C_Gamma(end, end) + (log_dM_R / (2*1i*pi));
end