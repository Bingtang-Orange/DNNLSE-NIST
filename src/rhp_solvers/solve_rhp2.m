function [q_val, info] = solve_rhp2(xev, t, q0, L_space, L_spectral, n_ray, n_real)
% SOLVE_RHP2 Solves deformed RHP II with dynamic scaling and zero-sum diagnostic.
%
% Syntax:
%   [q_val, info] = solve_rhp2(xev, t, q0)
%   [q_val, info] = solve_rhp2(xev, t, q0, L_space, L_spectral, n_ray, n_real)
%
% Inputs:
%   xev, t       - Spatial and temporal coordinates (x, t) with t > 0
%   q0           - Initial potential handle or Chebfun object
%   L_space      - Spatial truncation length for direct scattering (default: 20.0)
%   L_spectral   - Truncation radius of real spectral axis Gamma_2 (default: 8.0)
%   n_ray        - Collocation nodes per diagonal steepest descent ray (default: 24)
%   n_real       - Collocation nodes on truncated real segment Gamma_2 (default: 64)
%
% Outputs:
%   q_val        - Complex reconstructed potential value q(x, t)
%   info         - Diagnostic struct reporting saddle point lambda0, lens size d_ray,
%                  and the zero-sum compatibility residual max_k |sum +/- U_j(lambda0)|
%
% Description:
%   Implements the 5-contour steepest descent deformation of RHP II around the
%   stationary phase point lambda0 = x / (4*t). Dynamically scales the ray length
%   via d_ray(t) ~ O(t^{-1/2}), solves the global Cauchy singular system via
%   assemble_global, and reconstructs the potential via Clenshaw-Curtis integration.
%
% References:
%   Section 4.1 and Section 5.2.3.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    if t <= 1e-12
        error('solve_rhp2: Time t must be strictly positive (t > 0).');
    end
    if nargin < 4 || isempty(L_space),    L_space = 20.0;    end
    if nargin < 5 || isempty(L_spectral), L_spectral = 8.0;  end
    if nargin < 6 || isempty(n_ray),      n_ray = 24;        end
    if nargin < 7 || isempty(n_real),     n_real = 64;       end

    % 1. Stationary phase saddle point and dynamic steepest-descent ray scaling
    lambda0 = xev / (4 * t);
    d_ray = max(min(2.5, 2.5 / sqrt(t)), 0.05); % Matches Gaussian phase decay width

    v_SE = d_ray * exp(-1i * pi/4);
    v_NE = d_ray * exp( 1i * pi/4);
    v_NW = -conj(v_NE);
    v_SW = -conj(v_SE);

    end_real = max(lambda0 + 2.0, L_spectral);

    % Define 5-contour steepest descent deformation geometry
    defs = {
        create_curve_def('interval', lambda0, lambda0 + v_SE); ... % Gamma_1: L factor
        create_curve_def('interval', lambda0, end_real);       ... % Gamma_2: D factor
        create_curve_def('interval', lambda0, lambda0 + v_NE); ... % Gamma_3: U factor
        create_curve_def('interval', lambda0, lambda0 + v_NW); ... % Gamma_4: B^{-1} factor
        create_curve_def('interval', lambda0, lambda0 + v_SW)      % Gamma_5: A^{-1} factor
    };

    n_all = [n_ray, n_real, n_ray, n_ray, n_ray];
    N_total = sum(n_all);
    block_start = cumsum([0, n_all(1:end-1)]) + 1;
    block_end   = cumsum(n_all);

    % 2. Collocation mapping and direct scattering evaluation
    lambda_cell = cell(5, 1);
    x_cell = cell(5, 1);
    for k = 1:5
        xk = flipud(cos(pi * (0:n_all(k)-1)' / (n_all(k) - 1)));
        x_cell{k} = xk;
        pts = curve_mapping(defs{k}, xk);
        lambda_cell{k} = pts{1};
    end
    x_all = vertcat(x_cell{:});
    lambda_all = vertcat(lambda_cell{:});

    [rho_all, rho_hat_all] = compute_reflection(lambda_all, q0, L_space);

    % 3. Spatio-temporal jump matrix assembly
    G11 = ones(N_total, 1);  G12 = zeros(N_total, 1);
    G21 = zeros(N_total, 1); G22 = ones(N_total, 1);
    th_all = -lambda_all * xev + 2 * (lambda_all.^2) * t;

    for k = 1:5
        idx = block_start(k):block_end(k);
        th  = th_all(idx);
        r   = rho_all(idx);
        rh  = rho_hat_all(idx);
        det_fac = 1 - r .* rh;

        switch k
            case 1, G21(idx) = (rh .* exp(-2i * th)) ./ det_fac;
            case 2, G11(idx) = det_fac; G22(idx) = 1 ./ det_fac;
            case 3, G12(idx) = (-r .* exp(2i * th)) ./ det_fac;
            case 4, G21(idx) = -rh .* exp(-2i * th);
            case 5, G12(idx) =  r .* exp(2i * th);
        end
    end

    % 4. Solve global singular integral system
    [Cplus, Cminus] = assemble_global(defs, x_all, n_all);
    A11 = Cplus - diag(G11) * Cminus;
    A12 = -diag(G21) * Cminus;
    A21 = -diag(G12) * Cminus;
    A22 = Cplus - diag(G22) * Cminus;

    anss = [A11, A12; A21, A22] \ [G11 - 1; G12];
    U11 = anss(1:N_total);
    U12 = anss(N_total + 1 : 2*N_total);

    % 5. Zero-sum junction condition diagnostic at lambda0
    % All 5 contours originate at lambda0 (corresponding to node 1 of each block)
    idx_lambda0 = block_start;
    sum_U11 = sum(U11(idx_lambda0));
    sum_U12 = sum(U12(idx_lambda0));

    info.lambda0         = lambda0;
    info.d_ray           = d_ray;
    info.zero_sum_res_12 = abs(sum_U12);
    info.zero_sum_res_11 = abs(sum_U11);
    info.zero_sum_max    = max(abs(sum_U11), abs(sum_U12));

    % 6. Potential reconstruction via Clenshaw-Curtis quadrature
    q_val = 0;
    for k = 1:5
        idx = block_start(k):block_end(k);
        wk = clencurt(n_all(k));
        dz_dx = (defs{k}.b - defs{k}.a) / 2;
        q_val = q_val - (1 / pi) * (dz_dx * sum(wk(:) .* U12(idx)));
    end
end