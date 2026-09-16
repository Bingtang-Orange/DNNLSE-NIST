function S = compute_S(lambda, pot)
% COMPUTE_S Solves the spatial Lax pair ODEs to assemble the scattering matrix.
%
% Syntax:
%   S = compute_S(lambda, pot)
%
% Inputs:
%   lambda - Scalar spectral parameter
%   pot    - Struct containing precomputed subdomain potentials and domain length L
%
% Outputs:
%   S      - 2x2 complex scattering matrix S(lambda)
%
% References:
%   Defocusing nonlocal NLS Lax pair ODE collocation, Section 3.1, Eqs. (20)-(22).

    L    = pot.L;
    q_L  = pot.q_left;   qc_L = pot.q_neg_c_left;
    q_R  = pot.q_right;  qc_R = pot.q_neg_c_right;

    % ---- Left subdomain [-L, 0]: First-order ODE system for mu_1 ----
    % Column 1: [a1; c1] with boundary condition [a1(-L); c1(-L)] = [0; 0]
    L1a = chebop(-L, 0);
    L1a.op = @(x, a, c) [ diff(a) - q_L*c; ...
                          diff(c) - qc_L*a - 2i*lambda*c ];
    L1a.lbc = @(a, c) [a; c];
    U1a = L1a \ [0; qc_L];

    % Column 2: [b1; d1] with boundary condition [b1(-L); d1(-L)] = [0; 0]
    L1b = chebop(-L, 0);
    L1b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_L*d; ...
                          diff(d) - qc_L*b ];
    L1b.lbc = @(b, d) [b; d];
    U1b = L1b \ [q_L; 0];

    % ---- Right subdomain [0, L]: First-order ODE system for mu_2 ----
    % Column 1: [a2; c2] with boundary condition [a2(L); c2(L)] = [0; 0]
    L2a = chebop(0, L);
    L2a.op = @(x, a, c) [ diff(a) - q_R*c; ...
                          diff(c) - qc_R*a - 2i*lambda*c ];
    L2a.rbc = @(a, c) [a; c];
    U2a = L2a \ [0; qc_R];

    % Column 2: [b2; d2] with boundary condition [b2(L); d2(L)] = [0; 0]
    L2b = chebop(0, L);
    L2b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_R*d; ...
                          diff(d) - qc_R*b ];
    L2b.rbc = @(b, d) [b; d];
    U2b = L2b \ [q_R; 0];

    % ---- Extract solution values at the interface x = 0 ----
    mu1_0 = [U1a{1}(0), U1b{1}(0);
             U1a{2}(0), U1b{2}(0)];
    mu2_0 = [U2a{1}(0), U2b{1}(0);
             U2a{2}(0), U2b{2}(0)];

    % ---- Assemble scattering matrix S(lambda) via connection formula ----
    I = eye(2);
    S = (mu2_0 + I) \ (mu1_0 + I);
end