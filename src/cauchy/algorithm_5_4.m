function [mu_L, mu_R] = algorithm_5_4(n)
%ALGORITHM_5_4 Compute endpoint singular boundary correction vectors.
%   [MU_L, MU_R] = ALGORITHM_5_4(N) computes the discrete boundary limit
%   vectors for the Cauchy transform on the reference interval [-1, 1].
%
%   Inputs:
%       n    - Number of Chebyshev collocation nodes (scalar integer)
%
%   Outputs:
%       mu_L - Left endpoint correction row vector at x = -1 (size 1 x n)
%       mu_R - Right endpoint correction row vector at x = +1 (size 1 x n)

    mu_vals = zeros(1, n);
    for k = 1:n-1
        if mod(k, 2) ~= 0
            mu_vals(k+1) = mu_vals(k) + 1/k;
        else
            mu_vals(k+1) = mu_vals(k);
        end
    end

    mu_L = zeros(1, n);
    for k = 2:n
        mu_L(k) = ((-1)^k) * (mu_vals(k-1) + mu_vals(k));
    end
    mu_L = mu_L / (1i * pi);

    mu_R = zeros(1, n);
    for k = 2:n
        mu_R(k) = mu_vals(k-1) + mu_vals(k);
    end
    mu_R = mu_R / (1i * pi);
end