function [delta_val, chi_val] = compute_delta(lambda_targets, lambda0, q0, L_space, L_spectral, N_delta, which_side)
% COMPUTE_DELTA Solves the scalar RH problem for delta(lambda) on the branch cut [lambda0, +inf).
%
% Syntax:
%   [delta_val, chi_val] = compute_delta(lambda_targets, lambda0, q0)
%   [delta_val, chi_val] = compute_delta(lambda_targets, lambda0, q0, L_space, L_spectral, N_delta, which_side)
%
% Inputs:
%   lambda_targets - Complex evaluation coordinates (column vector)
%   lambda0        - Stationary phase saddle point x / (4*t)
%   q0             - Initial potential function handle or Chebfun object
%   L_space        - Spatial domain truncation half-length (default: 20.0)
%   L_spectral     - Spectral support half-width of the reflection coefficient (default: 8.0)
%   N_delta        - Number of Chebyshev nodes on branch cut [lambda0, L_spectral] (default: 64)
%   which_side     - Boundary limit orientation: '+' (upper), '-' (lower), or a char vector
%
% Outputs:
%   delta_val      - Evaluated scalar function delta(lambda)
%   chi_val        - Complex phase exponent chi(lambda) = log(delta(lambda))
%
% Description:
%   Solves the jump relation delta_+ = delta_- * (1 - rho*rho_hat) across [lambda0, +inf).
%   The semi-infinite cut is mapped to [-1, 1], expanded via Chebyshev polynomials,
%   and evaluated off the cut via the exterior Joukowsky conformal mapping.
%
% References:
%   Section 4.1 and Section 4.2.

    if nargin < 4 || isempty(L_space),    L_space = 20.0;    end
    if nargin < 5 || isempty(L_spectral), L_spectral = 8.0;  end
    if nargin < 6 || isempty(N_delta),    N_delta = 64;      end
    if nargin < 7 || isempty(which_side), which_side = '+';  end

    lambda_targets = lambda_targets(:);
    M_pts = length(lambda_targets);

    % 1. Rapid exit if saddle point lies beyond effective spectral support
    if lambda0 >= L_spectral - 1e-12
        chi_val = complex(zeros(M_pts, 1));
        delta_val = complex(ones(M_pts, 1));
        return;
    end

    % 2. Chebyshev polynomial expansion on the finite cut segment [lambda0, L_spectral]
    a = lambda0;
    b = max(lambda0 + 1.0, L_spectral);

    x_I = flipud(cos(pi * (0:N_delta-1)' / (N_delta - 1))); % Reference grid [-1, 1]
    xi_nodes = ((b - a) / 2) * x_I + ((b + a) / 2);         % Physical grid [a, b]

    % Compute reflection coefficients along the branch cut
    [rho_xi, rho_hat_xi] = compute_reflection(xi_nodes, q0, L_space);
    det_factor = 1.0 - rho_xi .* rho_hat_xi;
    det_factor(det_factor <= 1e-15) = 1e-15; % Regularize near-zero factors
    h_nodes = log(det_factor);

    % Chebyshev expansion coefficients
    F = inv(cheby(x_I, N_delta));
    c_coeffs = F * h_nodes;

    % 3. Map evaluation targets to reference conformal plane w
    w = (2 * lambda_targets - (b + a)) / (b - a);

    tol_zero = 1e-12;
    tol_edge = 1e-8;

    is_real_axis  = abs(imag(w)) < tol_zero;
    is_inside_cut = (real(w) > -1.0 + tol_edge) & (real(w) < 1.0 - tol_edge);
    on_cut_mask   = is_real_axis & is_inside_cut;

    % 4. Handle boundary limits for evaluation points residing on the cut
    if any(on_cut_mask)
        w_real = real(w(on_cut_mask));
        if numel(which_side) == M_pts
            side_cut = which_side(on_cut_mask);
            sgn = ones(sum(on_cut_mask), 1);
            sgn(side_cut == '-') = -1;
            w(on_cut_mask) = w_real + 1i * sgn * 1e-13;
        elseif strcmp(which_side, '-')
            w(on_cut_mask) = w_real - 1i * 1e-13;
        else
            w(on_cut_mask) = w_real + 1i * 1e-13;
        end
    end

    % Numerical boundary clipping near branch singularities
    w(real(w) <= -1.0 + tol_edge & is_real_axis) = -1.0 + tol_edge;
    w(real(w) >=  1.0 - tol_edge & is_real_axis) =  1.0 - tol_edge;

    % 5. Joukowsky inverse map and global Cauchy basis matrix evaluation
    z = w - sqrt(w - 1) .* sqrt(w + 1);
    Psi_z = psi_row(z, N_delta);

    chi_val   = Psi_z * c_coeffs;
    delta_val = exp(chi_val);
end