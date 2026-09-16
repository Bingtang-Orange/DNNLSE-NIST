function [mu_L, mu_R] = algorithm_5_4(n)
% ALGORITHM_5_4 Generates endpoint singular boundary correction vectors mu_L and mu_R.
%
% Syntax:
%   [mu_L, mu_R] = algorithm_5_4(n)
%
% Inputs:
%   n    - Polynomial basis dimension / number of Chebyshev collocation nodes
%
% Outputs:
%   mu_L - Left endpoint singular correction row vector at x = -1 (size 1 x n)
%   mu_R - Right endpoint singular correction row vector at x = +1 (size 1 x n)
%
% Description:
%   Implements Algorithm 5.4 of Olver (2012). Computes the discrete boundary limit
%   vectors for the Cauchy transform on the reference interval [-1, 1]. These vectors
%   are essential for rank-one regularization at contour junctions and rays to eliminate
%   discrete corner blow-up.
%
% References:
%   S. Olver, Numer. Math. 122 (2012), Algorithm 5.4.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % 1. Generate the intermediate scalar recurrence sequence mu_vals
    mu_vals = zeros(1, n);
    for k = 1:n-1
        if mod(k, 2) ~= 0
            mu_vals(k+1) = mu_vals(k) + 1/k;
        else
            mu_vals(k+1) = mu_vals(k);
        end
    end

    % 2. Assemble left boundary vector mu_L (corresponding to endpoint x = -1)
    mu_L = zeros(1, n);
    for k = 2:n
        mu_L(k) = ((-1)^k) * (mu_vals(k-1) + mu_vals(k));
    end
    mu_L = mu_L / (1i * pi);

    % 3. Assemble right boundary vector mu_R (corresponding to endpoint x = +1)
    mu_R = zeros(1, n);
    for k = 2:n
        mu_R(k) = mu_vals(k-1) + mu_vals(k);
    end
    mu_R = mu_R / (1i * pi);
end