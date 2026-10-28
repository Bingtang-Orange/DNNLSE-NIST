function C = cauchy_def_5_7(x, n, def_src, def_tgt)
%CAUCHY_DEF_5_7 Evaluate Cauchy operator for disjoint contours via conformal mapping.
%   C = CAUCHY_DEF_5_7(X, N, DEF_SRC, DEF_TGT) evaluates the discrete
%   Cauchy transform matrix between disjoint source and target contours.
%
%   Inputs:
%       x       - Collocation nodes on reference interval for target contour
%       n       - Number of Chebyshev basis polynomials on source contour
%       def_src - Source contour definition structure
%       def_tgt - Target contour definition structure
%
%   Outputs:
%       C       - Discrete Cauchy transform matrix of size length(x) x n

    z_omega = curve_mapping(def_tgt, x);
    z_omega = z_omega{1};

    w = forward_map(def_src, z_omega);
    def_identity = create_curve_def('interval', -1, 1);

    C = cauchy_def_5_3(w, n, def_identity);
end