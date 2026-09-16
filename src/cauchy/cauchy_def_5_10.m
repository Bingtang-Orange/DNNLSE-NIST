function C = cauchy_def_5_10(x, n, def_src, def_tgt)
% CAUCHY_DEF_5_10 Evaluates Cauchy operator C[Gamma, Omega] for arc source and disjoint target.
%
% Syntax:
%   C = cauchy_def_5_10(x, n, def_src, def_tgt)
%
% Inputs:
%   x       - Collocation nodes on reference interval for target contour Omega
%   n       - Polynomial basis dimension on source arc Gamma
%   def_src - Source arc definition struct (with finite M_Gamma(inf))
%   def_tgt - Disjoint target contour definition struct Omega
%
% Outputs:
%   C       - Discrete Cauchy transform matrix of size length(x) x n
%
% Description:
%   Subtracts the rank-one infinity correction term from the leading operator
%   cauchy_def_5_7 to account for the finite mapping of infinity by the source arc.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.10.

    % 1. Evaluate leading Cauchy operator term via Definition 5.7
    C_base = cauchy_def_5_7(x, n, def_src, def_tgt);

    % 2. Extract finite conformal image of infinity w_inf = M_Gamma(inf)
    [~, info] = curve_mapping(def_src, 0);
    w_inf = info(1).M_inf;

    % 3. Evaluate Joukowsky conformal inverse and basis row
    z_psi = w_inf - sqrt(w_inf - 1) * sqrt(w_inf + 1);
    Psi_z = psi_row(z_psi, n);

    % 4. Compute Chebyshev basis transformation matrix F
    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));

    % 5. Apply rank-one correction
    sub = ones(size(C_base, 1), 1) * (Psi_z * F);
    C = C_base - sub;
end