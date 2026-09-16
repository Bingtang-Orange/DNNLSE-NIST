function C = cauchy_def_5_14(x, n, def_src, def_tgt)
% CAUCHY_DEF_5_14 Evaluates Cauchy operator C[Gamma, Omega] for ray source and connected target.
%
% Syntax:
%   C = cauchy_def_5_14(x, n, def_src, def_tgt)
%
% Inputs:
%   x       - Collocation nodes on target curve Omega (length n_Omega)
%   n       - Number of collocation nodes on source ray Gamma
%   def_src - Source ray definition struct (M_inf = +/- 1)
%   def_tgt - Target contour definition struct meeting Gamma at its finite endpoint
%
% Outputs:
%   C       - Discrete Cauchy transform matrix of size n_Omega x n
%
% Description:
%   Regularizes the Cauchy integral from a ray to a contour intersecting at its
%   finite endpoint by evaluating cauchy_def_5_8 on the extended (n+1) grid,
%   slicing out the column corresponding to infinity, and applying the rank-one
%   boundary correction to guarantee uniform L^2 operator boundedness.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.14.

    m = n + 1; % Extended basis dimension including infinity

    % 1. Evaluate base connected-contour matrix on extended dimension m (n_Omega x m)
    C_58 = cauchy_def_5_8(x, m, def_src, def_tgt);

    % 2. Identify M_inf of source ray
    if strcmp(def_src.direction, 'forward')
        M_inf = 1;
    else
        M_inf = -1;
    end

    % 3. Retrieve boundary limit vector mu and extended Chebyshev transformation matrix
    [mu_L, mu_R] = algorithm_5_4(m);
    x_I = flipud(cos(pi * (0:m-1)' / (m - 1)));
    F = inv(cheby(x_I, m));

    % 4. Slice out column at infinity and apply rank-one correction
    if M_inf == -1
        % M_inf = -1: Slice out column 1 (infinity); retain columns 2:m
        C_sliced = C_58(:, 2:m);
        F_gamma  = F(:, 2:m);
        mu_vec   = mu_L;
    else
        % M_inf = +1: Slice out column m (infinity); retain columns 1:n
        C_sliced = C_58(:, 1:n);
        F_gamma  = F(:, 1:n);
        mu_vec   = mu_R;
    end

    n_Omega = size(C_sliced, 1);
    sub = ones(n_Omega, 1) * (mu_vec * F_gamma);
    C = C_sliced - sub;
end