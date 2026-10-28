function C = cauchy_def_5_10(x, n, def_src, def_tgt)
%CAUCHY_DEF_5_10 Evaluate Cauchy operator for arc source and disjoint target.
%   C = CAUCHY_DEF_5_10(X, N, DEF_SRC, DEF_TGT) evaluates the discrete
%   Cauchy transform matrix from a source arc to a disjoint target contour.
%
%   Inputs:
%       x       - Collocation nodes on reference interval for target contour
%       n       - Polynomial basis dimension on source arc
%       def_src - Source arc definition structure
%       def_tgt - Disjoint target contour definition structure
%
%   Outputs:
%       C       - Discrete Cauchy transform matrix of size length(x) x n

    C_base = cauchy_def_5_7(x, n, def_src, def_tgt);

    [~, info] = curve_mapping(def_src, 0);
    w_inf = info(1).M_inf;

    z_psi = w_inf - sqrt(w_inf - 1) * sqrt(w_inf + 1);
    Psi_z = psi_row(z_psi, n);

    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));

    sub = ones(size(C_base, 1), 1) * (Psi_z * F);
    C = C_base - sub;
end