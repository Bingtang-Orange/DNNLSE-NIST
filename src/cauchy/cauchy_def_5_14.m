function C = cauchy_def_5_14(x, n, def_src, def_tgt)
%CAUCHY_DEF_5_14 Evaluate Cauchy operator for ray source and connected target.
%   C = CAUCHY_DEF_5_14(X, N, DEF_SRC, DEF_TGT) evaluates the discrete
%   Cauchy transform matrix from an unbounded ray to an intersecting target contour.
%
%   Inputs:
%       x       - Collocation nodes on target curve (length n_Omega)
%       n       - Number of collocation nodes on source ray
%       def_src - Source ray definition structure
%       def_tgt - Target contour definition structure meeting ray at its endpoint
%
%   Outputs:
%       C       - Discrete Cauchy transform matrix of size length(x) x n

    m = n + 1;

    C_58 = cauchy_def_5_8(x, m, def_src, def_tgt);

    if strcmp(def_src.direction, 'forward')
        M_inf = 1;
    else
        M_inf = -1;
    end

    [mu_L, mu_R] = algorithm_5_4(m);
    x_I = flipud(cos(pi * (0:m-1)' / (m - 1)));
    F = inv(cheby(x_I, m));

    if M_inf == -1
        C_sliced = C_58(:, 2:m);
        F_gamma  = F(:, 2:m);
        mu_vec   = mu_L;
    else
        C_sliced = C_58(:, 1:n);
        F_gamma  = F(:, 1:n);
        mu_vec   = mu_R;
    end

    n_Omega = size(C_sliced, 1);
    sub = ones(n_Omega, 1) * (mu_vec * F_gamma);
    C = C_sliced - sub;
end