function C = cauchy_def_5_5(x, n, def)
% CAUCHY_DEF_5_5 Computes Cauchy matrix C[I, Gamma] for connected contours meeting at +/- 1.
%
% Syntax:
%   C = cauchy_def_5_5(x, n, def)
%
% Inputs:
%   x   - Reference collocation nodes for target curve Gamma (length n_Gamma)
%   n   - Number of polynomial basis functions on source interval I
%   def - Geometric definition struct of target curve Gamma
%
% Outputs:
%   C   - Discrete Cauchy matrix C[I, Gamma] with corner singularity regularization
%
% Description:
%   Regularizes the boundary endpoint singularity where Gamma touches I at +/- 1
%   via explicit tangent angle correction and discrete boundary vector mu.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.5 and Algorithm 5.4.

    tol = 1e-12;

    % 1. Identify junction endpoint topology
    if strcmp(def.type, 'interval')
        z_start = def.a; z_end = def.b;
    elseif strcmp(def.type, 'ray')
        if strcmp(def.direction, 'forward')
            z_start = def.a; z_end = inf;
        else
            z_start = inf; z_end = def.a;
        end
    elseif strcmp(def.type, 'arc')
        z_start = def.a + def.r * exp(1i * def.theta1);
        z_end   = def.a + def.r * exp(1i * def.theta2);
    end

    if abs(z_start - 1) < tol
        case_id = 1; % Curve originates from +1 (right end of I)
    elseif abs(z_start + 1) < tol
        case_id = 2; % Curve originates from -1 (left end of I)
    elseif abs(z_end - 1) < tol
        case_id = 3; % Curve terminates at +1 (right end of I)
    elseif abs(z_end + 1) < tol
        case_id = 4; % Curve terminates at -1 (left end of I)
    else
        error('cauchy_def_5_5: Curve endpoint must terminate at either +1 or -1.');
    end

    % 2. Tangent angle and singular boundary correction vector
    [dz_dx, ~, ~] = curve_derivative(def, x);
    k_vec = 0:n-1;
    r_R = ones(1, n) / (2*1i*pi);
    r_L = -(-1).^k_vec / (2*1i*pi);

    switch case_id
        case 1
            theta = angle(dz_dx(1));
            r_vec = r_R; sign_val = +1;
        case 2
            theta = angle(dz_dx(1));
            r_vec = r_L; sign_val = -1;
        case 3
            theta = angle(dz_dx(end)) - pi;
            r_vec = r_R; sign_val = +1;
        case 4
            theta = angle(dz_dx(end)) - pi;
            r_vec = r_L; sign_val = -1;
    end

    % 3. Evaluate boundary limit vector mu via Algorithm 5.4
    [mu_L, mu_R] = algorithm_5_4(n);
    if case_id == 2 || case_id == 4
        mu_vec = mu_L;
    else
        mu_vec = mu_R;
    end

    factor = 1i * angle(sign_val * exp(1i*theta)) - log(2);
    mu_row = mu_vec + factor * r_vec;

    % 4. Map nodes to physical coordinates and isolate non-touching interior points
    if strcmp(def.type, 'ray')
        if strcmp(def.direction, 'backward')
            x_full = [-1; x];
        else
            x_full = [x; 1];
        end
    else
        x_full = x;
    end
    pts = curve_mapping(def, x_full);
    z_full = pts{1};

    % Exclude intersecting junction node
    if case_id == 1 || case_id == 2
        z_remain = z_full(2:end);
    else
        z_remain = z_full(1:end-1);
    end

    % 5. Evaluate basis on strictly exterior interior nodes (|z_psi| < 1)
    z_psi = z_remain - sqrt(z_remain - 1) .* sqrt(1 + z_remain);
    Psi_z = psi_row(z_psi, n);

    % 6. Assemble regularized matrix and apply Chebyshev basis transformation
    if case_id == 1 || case_id == 2
        M = [mu_row; Psi_z];
    else
        M = [Psi_z; mu_row];
    end

    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));
    C = M * F;
end