function C_Gamma = cauchy_def_5_6(x, n, def)
% CAUCHY_DEF_5_6 Implements discrete Cauchy operator C^+[Gamma, Gamma] for bounded curves.
%
% Syntax:
%   C_Gamma = cauchy_def_5_6(x, n, def)
%
% Inputs:
%   x       - Collocation nodes on the reference unit interval
%   n       - Number of Chebyshev collocation nodes
%   def     - Bounded curve definition struct satisfying M(inf) = inf
%
% Outputs:
%   C_Gamma - Modified discrete Cauchy operator matrix on Gamma
%
% Description:
%   Applies logarithmic endpoint metric derivatives along the diagonal to preserve
%   operator boundedness across non-unit conformal mappings.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.6.

    % 1. Compute reference Cauchy matrix on the unit interval
    [C_II, ~, ~, ~, ~] = cauchy_kernel_interval(x, n);

    % 2. Compute logarithmic derivatives of the forward mapping at endpoints
    [~, log_dM_L, log_dM_R] = curve_derivative(def, x);

    % 3. Apply diagonal endpoint regularizations
    C_Gamma = C_II;
    C_Gamma(1, 1)     = C_Gamma(1, 1)     - (log_dM_L / (2*1i*pi));
    C_Gamma(end, end) = C_Gamma(end, end) + (log_dM_R / (2*1i*pi));
end