function C_Gamma = cauchy_def_5_12(x, n, def)
%CAUCHY_DEF_5_12 Discrete Cauchy operator C^+ on unbounded rays.
%   C_GAMMA = CAUCHY_DEF_5_12(X, N, DEF) computes the discrete Cauchy
%   operator matrix on an unbounded ray.
%
%   Inputs:
%       x       - Collocation nodes on reference interval (excluding infinity)
%       n       - Number of collocation nodes on the ray
%       def     - Ray definition structure
%
%   Outputs:
%       C_Gamma - Discrete Cauchy operator matrix on the ray of size n x n

    if ~isfield(def, 'M_inf') || isempty(def.M_inf)
        if strcmp(def.direction, 'forward')
            def.M_inf = 1;
        else
            def.M_inf = -1;
        end
    end
    M_inf = def.M_inf;

    m = n + 1;
    x_I = flipud(cos(pi*(0:m-1)'/(m-1)));

    [C_II, ~, ~, mu_L, mu_R] = cauchy_kernel_interval(x_I, m);

    F = inv(cheby(x_I, m));

    [~, log_dM_L, log_dM_R] = curve_derivative(def, x);

    if M_inf == -1
        C_Gamma = C_II(2:m, 2:m);
        F_Gamma = F(:, 2:m);
        mu_vec  = mu_L;
        C_Gamma(end, end) = C_Gamma(end, end) + (log_dM_R / (2*1i*pi));
    else
        C_Gamma = C_II(1:n, 1:n);
        F_Gamma = F(:, 1:n);
        mu_vec  = mu_R;
        C_Gamma(1, 1) = C_Gamma(1, 1) - (log_dM_L / (2*1i*pi));
    end

    sub = ones(n, 1) * (mu_vec * F_Gamma);
    C_Gamma = C_Gamma - sub;
end