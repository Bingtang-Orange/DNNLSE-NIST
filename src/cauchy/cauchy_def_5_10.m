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

    [~, info] = curve_mapping(def_src, 0);
    w_inf = info(1).M_inf;
    points = curve_mapping(def_tgt, x);
    z = points{1};
    c = def_src.r * exp(1i*(def_src.theta1+def_src.theta2)/2);
    numer = w_inf*(z-def_src.a-c);
    denom = z-def_src.a+c;
    large = abs(denom) < 0.25*abs(numer);
    mapped = zeros(size(z));
    inverse_w = denom(large)./numer(large);
    mapped(large) = inverse_w./(1+sqrt(1-inverse_w.^2));
    w = numer(~large)./denom(~large);
    mapped(~large) = 1./(w+sqrt(w-1).*sqrt(w+1));

    z_psi = w_inf - sqrt(w_inf - 1) * sqrt(w_inf + 1);
    Psi_z = psi_row(z_psi, n);

    x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
    F = inv(cheby(x_I, n));
    C_base = psi_row(mapped,n)*F;

    sub = ones(size(C_base, 1), 1) * (Psi_z * F);
    C = C_base - sub;
end