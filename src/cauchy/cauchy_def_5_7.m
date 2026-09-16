function C = cauchy_def_5_7(x, n, def_src, def_tgt)
% CAUCHY_DEF_5_7 Evaluates Cauchy operator C[Gamma, Omega] for disjoint contours via conformal mapping.
%
% Syntax:
%   C = cauchy_def_5_7(x, n, def_src, def_tgt)
%
% Inputs:
%   x       - Collocation nodes on reference interval for target contour Omega
%   n       - Number of Chebyshev basis polynomials on source contour Gamma
%   def_src - Source contour definition struct Gamma (bounded interval or circular arc)
%   def_tgt - Target contour definition struct Omega (disjoint from Gamma)
%
% Outputs:
%   C       - Discrete Cauchy transform matrix of size length(x) x n
%
% Description:
%   Evaluates C[Gamma, Omega] = C[I, M_Gamma(Omega)] by mapping target observation
%   nodes z_Omega into the source conformal coordinate w = M_Gamma(z_Omega),
%   and applying the exterior Joukowsky kernel via cauchy_def_5_3.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.7.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % 1. Map reference nodes to target contour coordinates in the complex plane
    z_omega = curve_mapping(def_tgt, x);
    z_omega = z_omega{1};

    % 2. Map observation points to source conformal plane: w = M_Gamma(z_omega)
    w = forward_map(def_src, z_omega);

    % 3. Set up canonical interval identity definition on [-1, 1]
    def_identity = create_curve_def('interval', -1, 1);

    % 4. Evaluate via reference exterior Cauchy operator cauchy_def_5_3
    C = cauchy_def_5_3(w, n, def_identity);
end