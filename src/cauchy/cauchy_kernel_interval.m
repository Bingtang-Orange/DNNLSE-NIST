function [Cplus, r_L, r_R, mu_L, mu_R] = cauchy_kernel_interval(x, n)
%CAUCHY_KERNEL_INTERVAL Discrete Cauchy operator and endpoint singular data on [-1, 1].
%   [CPLUS, R_L, R_R, MU_L, MU_R] = CAUCHY_KERNEL_INTERVAL(X, N) computes
%   the discrete Cauchy matrix C^+ on the reference interval [-1, 1] together
%   with endpoint regularization vectors.
%
%   Inputs:
%       x         - Chebyshev-Lobatto collocation nodes in [-1, 1] (column vector)
%       n         - Number of collocation nodes / polynomial basis order
%
%   Outputs:
%       Cplus     - Discrete Cauchy matrix C^+[I, I] of size n x n
%       r_L, r_R  - Endpoint rank-one singular vectors (size 1 x n)
%       mu_L, mu_R- Endpoint boundary limit vectors (size 1 x n)

    T = cheby(x, n);
    F = inv(T);

    P = zeros(n, n);
    for j = 1:ceil(n/2)
        denom  = 2*j - 1;
        offset = 2*j - 1;

        len = n - offset;
        if len <= 0
            break;
        end

        vals = (2 / denom) * ones(1, len);
        vals(1) = 1 / denom;
        P = P + diag(vals, offset);
    end

    D = zeros(n, n);
    D(1, 1) = log(2)/2 + 1i*pi/2;
    D(n, n) = -log(2)/2 + 1i*pi/2;

    x_bar = x(2:end-1);
    z_down = x_bar - 1i * sqrt(1 - x_bar) .* sqrt(1 + x_bar);
    arctanh_z = atanh(z_down);
    for k = 2:n-1
        D(k, k) = -2 * arctanh_z(k-1);
    end

    Cplus = (D + T * P * F) / (1i*pi);

    k = 0:n-1;
    r_L = -(-1).^k / (2*1i*pi);
    r_R = ones(1, n) / (2*1i*pi);

    mu_vals = zeros(1, n);
    for k = 1:n-1
        if mod(k, 2) == 1
            mu_vals(k+1) = mu_vals(k) + 1/k;
        else
            mu_vals(k+1) = mu_vals(k);
        end
    end

    mu_L = zeros(1, n);
    for k = 2:n
        mu_L(k) = ((-1)^k) * (mu_vals(k-1) + mu_vals(k));
    end
    mu_L = mu_L / (1i*pi);

    mu_R = zeros(1, n);
    for k = 2:n
        mu_R(k) = mu_vals(k-1) + mu_vals(k);
    end
    mu_R = mu_R / (1i*pi);
end