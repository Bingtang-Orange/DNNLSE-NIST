function [rho, rho_hat] = compute_reflection(lambda, q0, L)
% COMPUTE_REFLECTION Computes reflection coefficients for the DNNLSE Lax pair.
%
% Syntax:
%   [rho, rho_hat] = compute_reflection(lambda, q0)
%   [rho, rho_hat] = compute_reflection(lambda, q0, L)
%
% Inputs:
%   lambda  - Vector or scalar of spectral parameters (real or complex)
%   q0      - Initial potential function handle or Chebfun object
%   L       - Domain truncation half-length (default: 40.0)
%
% Outputs:
%   rho     - Reflection coefficient rho(lambda) = s_12 / s_22
%   rho_hat - Nonlocal reflection coefficient rho_hat(lambda) = s_21 / s_11
%
% Description:
%   Solves the spatial Lax pair ODE systems on (-L, 0] and [0, L) via 
%   Chebyshev spectral collocation (Chebfun), as described in Section 3.1.

    if nargin < 3 || isempty(L)
        L = 40.0; % Truncation length ensuring machine-precision decay for Schwartz data
    end

    lambda = lambda(:);
    N = numel(lambda);
    rho     = complex(zeros(N, 1));
    rho_hat = complex(zeros(N, 1));

    % Construct Chebfun representation with split domain [-L, 0, L]
    if ~isa(q0, 'chebfun')
        try
            q0 = chebfun(q0, [-L, 0, L]);
        catch
            q0 = chebfun(q0, [-L, L], 'splitting', 'on');
        end
    end

    % Pre-sample potentials on left and right subdomains to accelerate loop
    x_left  = chebfun('x', [-L, 0]);
    x_right = chebfun('x', [0, L]);

    potentials.q_left        = q0(x_left);
    potentials.q_neg_c_left  = conj(q0(-x_left));
    potentials.q_right       = q0(x_right);
    potentials.q_neg_c_right = conj(q0(-x_right));
    potentials.L             = L;

    % Evaluate scattering data across the spectral grid
    for k = 1:N
        S = compute_S(lambda(k), potentials);
        rho(k)     = S(1,2) / S(2,2);
        rho_hat(k) = S(2,1) / S(1,1);
    end
end