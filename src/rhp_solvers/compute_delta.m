function [delta_val, chi_val] = compute_delta(lambda_targets, lambda0, q0, L_space, L_spectral, N_delta, which_side)
%COMPUTE_DELTA Solve the scalar Riemann-Hilbert problem on [lambda0, +inf).
%   [DELTA_VAL, CHI_VAL] = COMPUTE_DELTA(LAMBDA_TARGETS, LAMBDA0, Q0)
%   [DELTA_VAL, CHI_VAL] = COMPUTE_DELTA(LAMBDA_TARGETS, LAMBDA0, Q0, L_SPACE, L_SPECTRAL, N_DELTA, WHICH_SIDE)
%   evaluates the scalar function DELTA_VAL and exponent CHI_VAL at points
%   LAMBDA_TARGETS across the branch cut [LAMBDA0, +inf).
%
%   Inputs:
%       lambda_targets - Complex evaluation coordinates (column vector)
%       lambda0        - Stationary phase saddle point x / (4*t)
%       q0             - Initial potential function handle or Chebfun object
%       L_space        - Spatial domain truncation half-length (optional, default: 20.0)
%       L_spectral     - Spectral domain truncation parameter (optional, default: 8.0)
%       N_delta        - Number of Chebyshev nodes on the cut (optional, default: 64)
%       which_side     - Boundary limit orientation: '+' or '-' (optional, default: '+')
%
%   Outputs:
%       delta_val      - Evaluated scalar function delta(lambda)
%       chi_val        - Complex phase exponent chi(lambda) = log(delta(lambda))

    if nargin < 4 || isempty(L_space),    L_space = 20.0;    end
    if nargin < 5 || isempty(L_spectral), L_spectral = 8.0;  end
    if nargin < 6 || isempty(N_delta),    N_delta = 64;      end
    if nargin < 7 || isempty(which_side), which_side = '+';  end

    lambda_targets = lambda_targets(:);
    M_pts = length(lambda_targets);

    if lambda0 >= L_spectral - 1e-12
        chi_val   = complex(zeros(M_pts, 1));
        delta_val = complex(ones(M_pts, 1));
        return;
    end

    a = lambda0;
    b = max(lambda0 + 1.0, L_spectral);

    x_I = flipud(cos(pi * (0:N_delta-1)' / (N_delta - 1)));
    xi_nodes = ((b - a) / 2) * x_I + ((b + a) / 2);

    [rho_xi, rho_hat_xi] = compute_reflection(xi_nodes, q0, L_space);
    det_factor = 1.0 - rho_xi .* rho_hat_xi;

    tol_singularity = 1e-14;
    if any(abs(det_factor) < tol_singularity)
        error('compute_delta: Singularity detected: 1 - rho*rho_hat vanishes on the branch cut.');
    end

    theta_continuous = unwrap(angle(det_factor));
    h_nodes = log(abs(det_factor)) + 1i * theta_continuous;

    F = inv(cheby(x_I, N_delta));
    c_coeffs = F * h_nodes;

    w = (2 * lambda_targets - (b + a)) / (b - a);

    tol_zero = 1e-12;
    tol_edge = 1e-8;

    is_real_axis  = abs(imag(w)) < tol_zero;
    is_inside_cut = (real(w) > -1.0 + tol_edge) & (real(w) < 1.0 - tol_edge);
    on_cut_mask   = is_real_axis & is_inside_cut;

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

    near_left_branch  = is_real_axis & (real(w) >= -1.0) & (real(w) <= -1.0 + tol_edge);
    near_right_branch = is_real_axis & (real(w) <=  1.0) & (real(w) >=  1.0 - tol_edge);
    w(near_left_branch)  = -1.0 + tol_edge;
    w(near_right_branch) =  1.0 - tol_edge;

    z = w - sqrt(w - 1) .* sqrt(w + 1);

    abs_z = abs(z);
    overflow_mask = (abs_z >= 1.0);
    if any(overflow_mask)
        z(overflow_mask) = z(overflow_mask) ./ (abs_z(overflow_mask) + eps) * (1.0 - 1e-14);
    end

    Psi_z = psi_row(z, N_delta);

    chi_val   = Psi_z * c_coeffs;
    delta_val = exp(chi_val);
end