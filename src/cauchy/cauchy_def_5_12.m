function C_Gamma = cauchy_def_5_12(x, n, def)
% CAUCHY_DEF_5_12 Discrete Cauchy operator C^+[Gamma, Gamma] for unbounded rays.
%
% Syntax:
%   C_Gamma = cauchy_def_5_12(x, n, def)
%
% Inputs:
%   x       - Collocation nodes on the reference interval (excluding infinity, length n)
%   n       - Number of collocation nodes on the ray
%   def     - Ray definition struct (satisfying M_inf = +/- 1)
%
% Outputs:
%   C_Gamma - Discrete Cauchy operator matrix on the ray of size n x n
%
% Description:
%   Constructs the Cauchy operator on an unbounded ray via an extended (n+1)-point
%   compactification on [-1, 1], removes the row/column corresponding to infinity,
%   applies the metric derivative correction at the finite endpoint, and regularizes
%   the operator via a rank-one projection using the singular boundary vector mu.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.12.

    % 1. Verify and assign the conformal image of infinity M_inf
    if ~isfield(def, 'M_inf') || isempty(def.M_inf)
        if strcmp(def.direction, 'forward')
            def.M_inf = 1;
        else
            def.M_inf = -1;
        end
    end
    M_inf = def.M_inf;

    % 2. Construct extended unit interval grid (dimension m = n + 1)
    m = n + 1;
    x_I = flipud(cos(pi*(0:m-1)'/(m-1))); % Ascending order on [-1, 1]

    % 3. Evaluate reference Cauchy matrix and boundary limit vectors on [-1, 1]
    [C_II, ~, ~, mu_L, mu_R] = cauchy_kernel_interval(x_I, m);

    % 4. Compute Chebyshev transformation matrix F on extended grid
    F = inv(cheby(x_I, m));

    % 5. Compute logarithmic metric derivatives at endpoints
    [~, log_dM_L, log_dM_R] = curve_derivative(def, x);

    % 6. Slice matrix and apply finite endpoint metric derivative correction
    if M_inf == -1
        % Infinity mapped to -1: slice out first row and column; keep right endpoint z_R
        C_Gamma = C_II(2:m, 2:m);
        F_Gamma = F(:, 2:m);
        mu_vec  = mu_L;

        % Apply metric correction to the bottom-right diagonal entry
        C_Gamma(end, end) = C_Gamma(end, end) + (log_dM_R / (2*1i*pi));
    else
        % Infinity mapped to +1: slice out last row and column; keep left endpoint z_L
        C_Gamma = C_II(1:n, 1:n);
        F_Gamma = F(:, 1:n);
        mu_vec  = mu_R;

        % Apply metric correction to the top-left diagonal entry
        C_Gamma(1, 1) = C_Gamma(1, 1) - (log_dM_L / (2*1i*pi));
    end

    % 7. Apply rank-one regularization: - ones(n, 1) * (mu_vec * F_Gamma)
    sub = ones(n, 1) * (mu_vec * F_Gamma);
    C_Gamma = C_Gamma - sub;
end