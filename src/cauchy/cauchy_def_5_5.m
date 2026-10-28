function C = cauchy_def_5_5(x, n, def)
%CAUCHY_DEF_5_5 Compute Cauchy matrix for connected contours meeting at +/- 1.
%   C = CAUCHY_DEF_5_5(X, N, DEF) computes the discrete Cauchy matrix
%   from the reference unit interval [-1, 1] to a connected target curve.
%
%   Inputs:
%       x   - Reference collocation nodes for target curve (length n_Gamma)
%       n   - Number of polynomial basis functions on source interval
%       def - Geometric definition structure of target curve
%
%   Outputs:
%       C   - Discrete Cauchy matrix of size length(x) x n

    tol = 1e-12;

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
        case_id = 1;
    elseif abs(z_start + 1) < tol
        case_id = 2;
    elseif abs(z_end - 1) < tol
        case_id = 3;
    elseif abs(z_end + 1) < tol
        case_id = 4;
    else
        error('cauchy_def_5_5: Curve endpoint must terminate at either +1 or -1.');
    end

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

    [mu_L, mu_R] = algorithm_5_4(n);
    if case_id == 2 || case_id == 4
        mu_vec = mu_L;
    else
        mu_vec = mu_R;
    end

    factor = 1i * angle(sign_val * exp(1i*theta)) - log(2);
    mu_row = mu_vec + factor * r_vec;

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

    if case_id == 1 || case_id == 2
        z_remain = z_full(2:end);
    else
        z_remain = z_full(1:end-1);
    end

    z_psi = z_remain - sqrt(z_remain - 1) .* sqrt(1 + z_remain);
    Psi_z = psi_row(z_psi, n);

    if case_id == 1 || case_id == 2
        M = [mu_row; Psi_z];
    else
        M = [Psi_z; mu_row];
    end

    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));
    C = M * F;
end