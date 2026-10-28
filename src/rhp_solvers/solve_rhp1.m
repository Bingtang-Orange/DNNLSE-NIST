function q_val = solve_rhp1(xev, t, lambda, rho, rho_hat, L, n)
%SOLVE_RHP1 Solve RHP I on the truncated real axis [-L, L].
%   Q_VAL = SOLVE_RHP1(XEV, T, LAMBDA, RHO, RHO_HAT, L, N) reconstructs
%   the complex potential Q_VAL at space-time point (XEV, T) by solving
%   the undeformed Riemann-Hilbert problem on [-L, L].
%
%   Inputs:
%       xev, t       - Spatial and temporal evaluation coordinates
%       lambda       - Collocation nodes on the truncated spectral axis [-L, L]
%       rho, rho_hat - Reflection coefficients evaluated on the spectral grid
%       L            - Spectral truncation half-width of the real axis
%       n            - Number of Chebyshev collocation nodes
%
%   Outputs:
%       q_val        - Complex reconstructed potential value q(x, t)

    x_I = flipud(cos(pi * (0:n-1)' / (n - 1)));
    def = create_curve_def('interval', -L, L);

    Cplus  = cauchy_def_5_6(x_I, n, def);
    Cminus = Cplus - eye(n);

    th  = -lambda * xev + 2 * (lambda.^2) * t;
    T11 = 1 - rho .* rho_hat;
    T12 = -rho .* exp(2i * th);
    T21 =  rho_hat .* exp(-2i * th);
    T22 = ones(n, 1);

    A11 = Cplus - diag(T11) * Cminus;
    A12 = -diag(T21) * Cminus;
    A21 = -diag(T12) * Cminus;
    A22 = Cplus - diag(T22) * Cminus;

    Coeffs = [ A11, A12 ;
               A21, A22 ];

    rhs = [ T11 - 1 ; 
            T12     ];

    sol = Coeffs \ rhs;
    U12 = sol(n+1 : 2*n);

    w = clencurt(n);
    q_val = -(L / pi) * sum(w(:) .* U12(:));
end