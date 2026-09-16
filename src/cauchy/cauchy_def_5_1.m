function C = cauchy_def_5_1(x, n)
% CAUCHY_DEF_5_1 Implements the discrete Cauchy operator C^+[I, I] on the unit interval.
%
% Syntax:
%   C = cauchy_def_5_1(x, n)
%
% Inputs:
%   x - Collocation nodes on the unit interval [-1, 1] (including endpoints)
%   n - Number of Chebyshev collocation nodes
%
% Outputs:
%   C - Discrete Cauchy matrix C^+[I, I] of size n x n
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.1.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % Direct evaluation of C^+ on the reference interval [-1, 1]
    [C, ~, ~, ~, ~] = cauchy_kernel_interval(x, n);
end