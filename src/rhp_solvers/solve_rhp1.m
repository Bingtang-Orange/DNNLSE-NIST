function q_val = solve_rhp1(xev, t, lambda, rho, rho_hat, L, n)
% SOLVE_RHP1 Solves RHP I on the truncated real axis [-L, L] for short-time potential recovery.
%
% Syntax:
%   q_val = solve_rhp1(xev, t, lambda, rho, rho_hat, L, n)
%
% Inputs:
%   xev, t       - Spatial and temporal evaluation coordinates (x, t)
%   lambda       - Collocation nodes on the truncated spectral axis [-L, L] (length n)
%   rho, rho_hat - Reflection coefficients evaluated on the spectral grid
%   L            - Spectral truncation half-width of the real axis contour
%   n            - Number of Chebyshev collocation nodes
%
% Outputs:
%   q_val        - Complex reconstructed potential value q(x, t)
%
% Description:
%   Implements the discrete collocation solver for RHP I on the bounded segment
%   [-L, L] via Olver (2012) Definition 5.6. Solves the singular integral system
%   for the first-row discrete potential (U11, U12) and reconstructs q(x, t) via
%   Clenshaw-Curtis quadrature.
%
% References:
%   S. Olver, Numer. Math. 122 (2012), Section 5.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % 1. Construct mapped discrete Cauchy operators on [-L, L] (Definition 5.6)
    x_I = flipud(cos(pi * (0:n-1)' / (n - 1))); % Reference Chebyshev nodes on [-1, 1]
    def = create_curve_def('interval', -L, L);

    Cplus  = cauchy_def_5_6(x_I, n, def);
    Cminus = Cplus - eye(n);

    % 2. Assemble spatio-temporal jump matrix elements
    th  = -lambda * xev + 2 * (lambda.^2) * t;
    T11 = 1 - rho .* rho_hat;
    T12 = -rho .* exp(2i * th);
    T21 =  rho_hat .* exp(-2i * th);
    T22 = ones(n, 1);

    % 3. Assemble discrete block linear system for the first row: [U11; U12]
    A11 = Cplus - diag(T11) * Cminus;
    A12 = -diag(T21) * Cminus;
    A21 = -diag(T12) * Cminus;
    A22 = Cplus - diag(T22) * Cminus;

    Coeffs = [ A11, A12 ;
               A21, A22 ];

    rhs = [ T11 - 1 ; 
            T12     ];

    % Direct linear solve
    sol = Coeffs \ rhs;
    U12 = sol(n+1 : 2*n);

    % 4. Potential reconstruction via Clenshaw-Curtis quadrature
    % q(x, t) = 2i * lim_{lambda -> inf} (lambda * M)_{12} = -(L / pi) * int_{-1}^1 U12 dx
    w = clencurt(n);
    q_val = -(L / pi) * sum(w(:) .* U12(:));
end