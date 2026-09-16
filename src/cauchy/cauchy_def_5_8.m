function C = cauchy_def_5_8(x, n, def_src, def_tgt)
% CAUCHY_DEF_5_8 Implements discrete Cauchy operator C[Gamma, Omega] for connected contours.
%
% Syntax:
%   C = cauchy_def_5_8(x, n, def_src, def_tgt)
%
% Inputs:
%   x       - Collocation nodes on target contour Omega (length n_Omega)
%   n       - Number of polynomial basis functions on source contour Gamma
%   def_src - Source contour definition struct Gamma (bounded interval or arc)
%   def_tgt - Target contour definition struct Omega (meeting Gamma at an endpoint)
%
% Outputs:
%   C       - Discrete Cauchy transform matrix of size n_Omega x n
%
% Description:
%   Resolves endpoint singularity where Omega intersects Gamma via:
%   (i)   conformal chain-rule derivative calculations dw/dx = M_src'(z) * dz/dx;
%   (ii)  logarithmic phase corrections and singular boundary vectors mu;
%   (iii) boundary derivative rank-one corrections to preserve L^2 operator boundedness.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.8 and Algorithm 5.4.

    tol = 1e-12;
    n_Omega = length(x);

    % 1. Determine junction topology (Cases 1-4)
    [tgt_start, tgt_end] = get_curve_endpoints(def_tgt);

    w_start = forward_map(def_src, tgt_start);
    w_end   = forward_map(def_src, tgt_end);

    if abs(w_start - 1) < tol
        case_id = 1;  % Target starts at w = +1 (right end of source)
    elseif abs(w_start + 1) < tol
        case_id = 2;  % Target starts at w = -1 (left end of source)
    elseif abs(w_end - 1) < tol
        case_id = 3;  % Target terminates at w = +1 (right end of source)
    elseif abs(w_end + 1) < tol
        case_id = 4;  % Target terminates at w = -1 (left end of source)
    else
        error('cauchy_def_5_8: Target endpoint must connect to +/-1 in the conformal w-plane.');
    end

    % 2. Compute composite chain-rule derivatives dw/dx = M_src'(z) * dz/dx
    if strcmp(def_tgt.type, 'ray')
        if isfield(def_tgt, 'direction') && strcmp(def_tgt.direction, 'backward')
            x_tgt_full = [-1; x(:)];
        else
            x_tgt_full = [x(:); 1];
        end
    else
        x_tgt_full = x(:);
    end

    [dz_dx_tgt, ~, ~] = curve_derivative(def_tgt, x_tgt_full);
    [~, M_prime_start] = forward_map(def_src, tgt_start);
    [~, M_prime_end]   = forward_map(def_src, tgt_end);

    dw_dx_L = M_prime_start * dz_dx_tgt(1);
    dw_dx_R = M_prime_end   * dz_dx_tgt(end);

    % 3. Set tangent angles, correction vectors, and directional signs
    k_vec = 0:n-1;
    r_R = ones(1, n) / (2*1i*pi);
    r_L = -(-1).^k_vec / (2*1i*pi);

    switch case_id
        case 1
            theta = angle(dw_dx_L);
            r_vec = r_R; sign_val = +1;
        case 2
            theta = angle(dw_dx_L);
            r_vec = r_L; sign_val = -1;
        case 3
            theta = angle(dw_dx_R) - pi;
            r_vec = r_R; sign_val = +1;
        case 4
            theta = angle(dw_dx_R) - pi;
            r_vec = r_L; sign_val = -1;
    end

    % 4. Retrieve boundary limit vector mu via Algorithm 5.4
    [mu_L, mu_R] = get_mu_vectors(n);
    if case_id == 2 || case_id == 4
        mu_vec = mu_L;
    else
        mu_vec = mu_R;
    end

    % 5. Assemble regularized endpoint boundary row
    factor = 1i * angle(sign_val * exp(1i*theta)) - log(2);
    mu_row = mu_vec + factor * r_vec;

    % 6. Evaluate basis on strictly exterior interior nodes (|z_psi| < 1)
    pts_tgt = curve_mapping(def_tgt, x_tgt_full);
    z_full = pts_tgt{1};

    if case_id == 1 || case_id == 2
        z_remain = z_full(2:end);
    else
        z_remain = z_full(1:end-1);
    end

    w_remain = forward_map(def_src, z_remain);
    z_psi = w_remain - sqrt(w_remain - 1) .* sqrt(1 + w_remain);
    Psi_z = psi_row(z_psi, n);

    % 7. Assemble regularized matrix M
    if case_id == 1 || case_id == 2
        M = [mu_row; Psi_z];
    else
        M = [Psi_z; mu_row];
    end

    % 8. Apply inverse discrete Chebyshev transformation
    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));
    C = M * F;

    % 9. Apply endpoint metric derivative corrections (Definition 5.8)
    x_dummy = [-1; 1];
    [~, log_dM_L, log_dM_R] = curve_derivative(def_src, x_dummy);

    switch case_id
        case 1
            C(1, n) = C(1, n) + log_dM_R / (2*1i*pi);
        case 2
            C(1, 1) = C(1, 1) - log_dM_L / (2*1i*pi);
        case 3
            C(n_Omega, n) = C(n_Omega, n) + log_dM_R / (2*1i*pi);
        case 4
            C(n_Omega, 1) = C(n_Omega, 1) - log_dM_L / (2*1i*pi);
    end
end

% =========================================================================
% Helper Function 1: Extract Physical Curve Endpoints
% =========================================================================
function [z_start, z_end] = get_curve_endpoints(def)
    switch def.type
        case 'interval'
            z_start = def.a;
            z_end   = def.b;
        case 'ray'
            if isfield(def, 'direction') && strcmp(def.direction, 'backward')
                z_start = Inf;
                z_end   = def.a;
            else
                z_start = def.a;
                z_end   = Inf;
            end
        case 'arc'
            z_start = def.a + def.r * exp(1i * def.theta1);
            z_end   = def.a + def.r * exp(1i * def.theta2);
        otherwise
            error('Unknown curve type: %s', def.type);
    end
end

% =========================================================================
% Helper Function 2: Generate Singular Boundary Vectors mu (Algorithm 5.4)
% =========================================================================
function [mu_L, mu_R] = get_mu_vectors(n)
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
    mu_L = mu_L / (1i*pi);

    mu_R = zeros(1, n);
    for k = 2:n
        mu_R(k) = (mu_vals(k-1) + mu_vals(k));
    end
    mu_R = mu_R / (1i*pi);
end