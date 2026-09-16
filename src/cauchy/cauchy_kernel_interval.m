function [Cplus, r_L, r_R, mu_L, mu_R] = cauchy_kernel_interval(x, n)
% CAUCHY_KERNEL_INTERVAL Discrete Cauchy operator C^+ on [-1,1] and endpoint singular data.
%
% Syntax:
%   [Cplus, r_L, r_R, mu_L, mu_R] = cauchy_kernel_interval(x, n)
%
% Inputs:
%   x       - Chebyshev-Lobatto collocation nodes in [-1, 1] (column vector, length n)
%   n       - Number of collocation nodes / polynomial basis order
%
% Outputs:
%   Cplus   - Discrete Cauchy matrix C^+[I, I] of size n x n
%   r_L,r_R - Endpoint singular rank-one vectors r_L and r_R (size 1 x n)
%   mu_L,mu_R - Endpoint singular limit vectors mu_L and mu_R via Algorithm 5.4 (size 1 x n)
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.1 and Algorithm 5.4.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % 1. Evaluate Chebyshev polynomial matrix T and its inverse transformation F
    T = cheby(x, n);
    F = inv(T);

    % 2. Assemble strictly upper-triangular matrix P (Definition 5.1)
    P = zeros(n, n);
    for j = 1:ceil(n/2)
        denom  = 2*j - 1;
        offset = 2*j - 1;

        len = n - offset;
        if len <= 0
            break;
        end

        vals = (2 / denom) * ones(1, len);
        vals(1) = 1 / denom; % T_0 coefficient scaling
        P = P + diag(vals, offset);
    end

    % 3. Assemble diagonal evaluation matrix D
    D = zeros(n, n);
    
    % Left endpoint (x = -1)
    D(1, 1) = log(2)/2 + 1i*pi/2;
    
    % Right endpoint (x = 1)
    D(n, n) = -log(2)/2 + 1i*pi/2;
    
    % Interior collocation points
    x_bar = x(2:end-1);
    z_down = x_bar - 1i * sqrt(1 - x_bar) .* sqrt(1 + x_bar);
    arctanh_z = atanh(z_down);
    for k = 2:n-1
        D(k, k) = -2 * arctanh_z(k-1);
    end

    % 4. Assemble discrete Cauchy matrix C^+[I, I]
    Cplus = (D + T * P * F) / (1i*pi);

    % 5. Evaluate endpoint singular vectors r_L and r_R
    k = 0:n-1;
    r_L = -(-1).^k / (2*1i*pi);
    r_R = ones(1, n) / (2*1i*pi);

    % 6. Generate scalar sequence for boundary limits (Algorithm 5.4)
    mu_vals = zeros(1, n);
    for k = 1:n-1
        if mod(k, 2) == 1
            mu_vals(k+1) = mu_vals(k) + 1/k;
        else
            mu_vals(k+1) = mu_vals(k);
        end
    end

    % 7. Assemble left boundary vector mu_L (x = -1)
    mu_L = zeros(1, n);
    for k = 2:n
        mu_L(k) = ((-1)^k) * (mu_vals(k-1) + mu_vals(k));
    end
    mu_L = mu_L / (1i*pi);

    % 8. Assemble right boundary vector mu_R (x = +1)
    mu_R = zeros(1, n);
    for k = 2:n
        mu_R(k) = mu_vals(k-1) + mu_vals(k);
    end
    mu_R = mu_R / (1i*pi);
end