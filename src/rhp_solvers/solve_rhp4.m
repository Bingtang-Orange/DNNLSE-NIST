function [q_val, info] = solve_rhp4(xev, t, q0, L_space, L_spectral, n_box, n_ray, r_box_in, d_ray_in)
% SOLVE_RHP4 Solves the fully regularized 9-contour RHP IV with dynamic lens scaling.
%
% Syntax:
%   [q_val, info] = solve_rhp4(xev, t, q0)
%   [q_val, info] = solve_rhp4(xev, t, q0, L_space, L_spectral, n_box, n_ray, r_box_in, d_ray_in)
%
% Inputs:
%   xev, t       - Spatial and temporal evaluation coordinates (x, t) with t > 0
%   q0           - Initial potential handle or Chebfun object
%   L_space      - Spatial domain truncation length for direct scattering (default: 20.0)
%   L_spectral   - Spectral domain truncation for scalar delta evaluation (default: 8.0)
%   n_box        - Number of collocation points per box edge Sigma_1 - Sigma_5 (default: 24)
%   n_ray        - Number of collocation points per exterior ray Sigma_6 - Sigma_9 (default: 20)
%   r_box_in     - (Optional) User-specified box half-width
%   d_ray_in     - (Optional) User-specified ray truncation length
%
% Outputs:
%   q_val        - Complex reconstructed potential value q(x, t)
%   info         - Diagnostic struct containing geometric parameters and the maximum
%                  zero-sum junction residual across the 5 intersection vertices
%
% References:
%   Section 4.2 (RHP IV formulation), Section 4.3 (Scaling and zero-sum compatibility),
%   and Section 5.2.3 (Spectral convergence).

    if nargin < 4 || isempty(L_space),    L_space = 20.0;    end
    if nargin < 5 || isempty(L_spectral), L_spectral = 8.0;  end
    if nargin < 6 || isempty(n_box),      n_box = 24;        end
    if nargin < 7 || isempty(n_ray),      n_ray = 20;        end

    % 1. Stationary phase saddle point and dynamic geometric scaling
    lambda0 = xev / (4 * t);
    
    if nargin >= 8 && ~isempty(r_box_in)
        r_box = r_box_in;
    else
        % Equation (32): bounds Re(2i*theta) <= 1.0 to prevent jump matrix overflow
        r_box = max(min(0.25, 0.35 / sqrt(t)), 0.005);
    end
    
    if nargin >= 9 && ~isempty(d_ray_in)
        d_ray = d_ray_in;
    else
        % Equation (33): guarantees outer jump decay to machine precision
        d_ray = max(min(2.50, 2.50 / sqrt(t)), 0.020);
    end

    % Define local box vertices and real axis entry point
    P_E = lambda0 + r_box;
    z1  = lambda0 + r_box + 1i * r_box;
    z2  = lambda0 - r_box + 1i * r_box;
    z3  = lambda0 - r_box - 1i * r_box;
    z4  = lambda0 + r_box - 1i * r_box;

    % 9-contour local lens configuration (Sigma_1 to Sigma_9)
    defs = {
        create_curve_def('interval', P_E, z1); ...                             % Sigma_1
        create_curve_def('interval', z1, z2); ...                              % Sigma_2
        create_curve_def('interval', z2, z3); ...                              % Sigma_3
        create_curve_def('interval', z3, z4); ...                              % Sigma_4
        create_curve_def('interval', z4, P_E); ...                             % Sigma_5
        create_curve_def('interval', z4, z4 + d_ray*exp(-1i*pi/4)); ...        % Sigma_6 (SE)
        create_curve_def('interval', z1, z1 + d_ray*exp( 1i*pi/4)); ...        % Sigma_7 (NE)
        create_curve_def('interval', z2, z2 - conj(d_ray*exp( 1i*pi/4))); ... % Sigma_8 (NW)
        create_curve_def('interval', z3, z3 - conj(d_ray*exp(-1i*pi/4)))      % Sigma_9 (SW)
    };

    num_curves = 9;
    n_all = [repmat(n_box, 1, 5), repmat(n_ray, 1, 4)];
    N_total = sum(n_all);
    block_start = cumsum([0, n_all(1:end-1)]) + 1;
    block_end   = cumsum(n_all);

    % 2. Map nodes to physical complex contours
    lambda_cell = cell(num_curves, 1);
    x_cell = cell(num_curves, 1);
    for k = 1:num_curves
        nk = n_all(k);
        xk = flipud(cos(pi * (0:nk-1)' / (nk - 1)));
        x_cell{k} = xk;
        pts = curve_mapping(defs{k}, xk);
        lambda_cell{k} = pts{1};
    end
    x_all = vertcat(x_cell{:});
    lambda_all = vertcat(lambda_cell{:});

    % 3. Evaluate direct scattering data and regularizing scalar delta
    [rho_all, rho_hat_all] = compute_reflection(lambda_all, q0, L_space);

    % Vectorized boundary side tagging for branch cut crossings
    side_all = repmat('+', N_total, 1);
    for k = [4, 5, 6, 9]
        side_all(block_start(k):block_end(k)) = '-';
    end
    delta_all = compute_delta(lambda_all, lambda0, q0, L_space, L_spectral, 64, side_all);

    % 4. Assemble global jump matrix entries
    G11 = ones(N_total, 1);  G12 = zeros(N_total, 1);
    G21 = zeros(N_total, 1); G22 = ones(N_total, 1);
    th_all = -lambda_all * xev + 2 * (lambda_all.^2) * t;

    for k = 1:num_curves
        idx = block_start(k):block_end(k);
        th  = th_all(idx);
        r   = rho_all(idx);
        rh  = rho_hat_all(idx);
        dl  = delta_all(idx);
        det_fac = 1.0 - r .* rh;

        switch k
            case 1
                G11(idx) = dl ./ det_fac;
                G22(idx) = det_fac ./ dl;
            case 2
                G11(idx) = dl ./ det_fac;
                G12(idx) = dl .* r .* exp(2i * th);
                G22(idx) = det_fac ./ dl;
            case 3
                G11(idx) = dl ./ det_fac;
                G12(idx) = dl .* r .* exp(2i * th);
                G21(idx) = (rh .* exp(-2i * th)) ./ (dl .* det_fac);
                G22(idx) = 1.0 ./ dl;
            case 4
                G11(idx) = dl;
                G21(idx) = (rh .* exp(-2i * th)) ./ (dl .* det_fac);
                G22(idx) = 1.0 ./ dl;
            case 5
                G11(idx) = dl;
                G22(idx) = 1.0 ./ dl;
            case 6
                G21(idx) = (rh .* exp(-2i * th)) ./ ((dl.^2) .* det_fac);
            case 7
                G12(idx) = -(dl.^2) .* r .* exp(2i * th) ./ det_fac;
            case 8
                G21(idx) = -(rh .* exp(-2i * th)) ./ (dl.^2);
            case 9
                G12(idx) =  (dl.^2) .* r .* exp(2i * th);
        end
    end

    % 5. Solve the global discrete Cauchy singular system
    [Cplus, Cminus] = assemble_global(defs, x_all, n_all);
    A11 = Cplus - diag(G11) * Cminus;
    A12 = -diag(G21) * Cminus;
    A21 = -diag(G12) * Cminus;
    A22 = Cplus - diag(G22) * Cminus;

    anss = [A11, A12; A21, A22] \ [G11 - 1; G12];
    U11 = anss(1:N_total);
    U12 = anss(N_total + 1 : 2*N_total);

    % 6. Zero-sum condition diagnostic across the 5 intersection nodes
    U_start = @(c_idx) [U11(block_start(c_idx)); U12(block_start(c_idx))];
    U_end   = @(c_idx) [U11(block_end(c_idx));   U12(block_end(c_idx))];

    res_PE = norm(U_start(1) - U_end(5), inf);
    res_z1 = norm(-U_end(1) + U_start(2) + U_start(7), inf);
    res_z2 = norm(-U_end(2) + U_start(3) + U_start(8), inf);
    res_z3 = norm(-U_end(3) + U_start(4) + U_start(9), inf);
    res_z4 = norm(-U_end(4) + U_start(5) + U_start(6), inf);

    info.lambda0      = lambda0;
    info.r_box        = r_box;
    info.d_ray        = d_ray;
    info.zero_sum_max = max([res_PE, res_z1, res_z2, res_z3, res_z4]);

    % 7. Potential reconstruction via Clenshaw-Curtis integration
    q_val = 0;
    for k = 1:num_curves
        idx = block_start(k):block_end(k);
        wk = clencurt(n_all(k));
        dz_dx = (defs{k}.b - defs{k}.a) / 2;
        q_val = q_val - (1 / pi) * (dz_dx * sum(wk(:) .* U12(idx)));
    end
end