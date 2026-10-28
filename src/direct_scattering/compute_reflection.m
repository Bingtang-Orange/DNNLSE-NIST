function [rho, rho_hat] = compute_reflection(lambda, q0, L)
%COMPUTE_REFLECTION Compute reflection coefficients for the defocusing nonlocal NLSE.
%   [RHO, RHO_HAT] = COMPUTE_REFLECTION(LAMBDA, Q0)
%   [RHO, RHO_HAT] = COMPUTE_REFLECTION(LAMBDA, Q0, L) computes the reflection
%   coefficients RHO and RHO_HAT for the spatial Lax pair of the defocusing nonlocal
%   nonlinear Schrodinger equation using Chebyshev spectral collocation.
%
%   Inputs:
%       lambda  - Vector or scalar of spectral parameters (real or complex)
%       q0      - Initial potential function handle or Chebfun object
%       L       - Domain truncation half-length (optional, default: 40.0)
%
%   Outputs:
%       rho     - Reflection coefficient rho(lambda) = s_12 / s_22
%       rho_hat - Nonlocal reflection coefficient rho_hat(lambda) = s_21 / s_11

    if nargin < 3 || isempty(L)
        L = 40.0;
    end

    lambda = lambda(:);
    N = numel(lambda);
    rho     = complex(zeros(N, 1));
    rho_hat = complex(zeros(N, 1));

    if ~isa(q0, 'chebfun')
        try
            q0 = chebfun(q0, [-L, 0, L]);
        catch
            q0 = chebfun(q0, [-L, L], 'splitting', 'on');
        end
    end

    x_left  = chebfun('x', [-L, 0]);
    x_right = chebfun('x', [0, L]);

    potentials.q_left        = q0(x_left);
    potentials.q_neg_c_left  = conj(q0(-x_left));
    potentials.q_right       = q0(x_right);
    potentials.q_neg_c_right = conj(q0(-x_right));
    potentials.L             = L;

    for k = 1:N
        S = compute_S(lambda(k), potentials);
        rho(k)     = S(1,2) / S(2,2);
        rho_hat(k) = S(2,1) / S(1,1);
    end
end