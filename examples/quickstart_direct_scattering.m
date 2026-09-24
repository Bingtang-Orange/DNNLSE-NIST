% =========================================================================
% Script: direct_scattering_general.m
% Description: Standalone numerical direct scattering solver for generic
%              smooth potentials q(x) under the Defocusing Nonlocal NLS
%              (DNNLSE) framework via Chebyshev spectral collocation.
%
% Mathematical Formulation:
%   Solves the decoupled first-order ODE systems on [-L, 0] and [0, L]
%   to reconstruct the left and right Jost solutions at the origin:
%     S(lambda) = [phi_2(0, 0, lambda)]^{-1} * phi_1(0, 0, lambda)
%
% Reference:
%   "Numerical Inverse Scattering Transform for the Defocusing Nonlocal
%    Nonlinear Schrodinger Equation", Section 3.1.
% Repository: https://github.com/Bingtang-Orange/DNNLSE-NIST
% =========================================================================

clear; clc;
warning('off', 'chebfun:constructor');

%% 1. Computational Setup and Parameters
L = 20;                             % Domain truncation half-width
lambda = 0.5;                       % Spectral parameter (real or complex)

% Define potential profile q(x) on [-L, L]
% (Can be replaced with any smooth localized complex profile)
q_full = chebfun(@(x) exp(-x.^2), [-L, L], 'splitting', 'on');

%% 2. Precompute Piecewise Potentials on Half-Lines
% Pre-evaluating q(x) and its nonlocal reflection conj(q(-x)) on [-L, 0]
% and [0, L] significantly accelerates the Chebyshev collocation solver.
q_L  = chebfun(@(x) q_full(x),        [-L, 0]);
qc_L = chebfun(@(x) conj(q_full(-x)), [-L, 0]);

q_R  = chebfun(@(x) q_full(x),        [0, L]);
qc_R = chebfun(@(x) conj(q_full(-x)), [0, L]);

%% 3. Solve for Left Jost Solution on [-L, 0]
% Column 1: [a1; c1] satisfies system (19a) with mu1(-L) = 0
L1a = chebop(-L, 0);
L1a.op = @(x, a, c) [ diff(a) - q_L*c; ...
                      diff(c) - qc_L*a - 2i*lambda*c - qc_L ];
L1a.lbc = @(a, c) [a; c];
U1a = L1a \ 0;

% Column 2: [b1; d1] satisfies system (19b) with mu1(-L) = 0
L1b = chebop(-L, 0);
L1b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_L*d - q_L; ...
                      diff(d) - qc_L*b ];
L1b.lbc = @(b, d) [b; d];
U1b = L1b \ 0;

%% 4. Solve for Right Jost Solution on [0, L]
% Column 1: [a2; c2] satisfies system (19a) with mu2(L) = 0
L2a = chebop(0, L);
L2a.op = @(x, a, c) [ diff(a) - q_R*c; ...
                      diff(c) - qc_R*a - 2i*lambda*c - qc_R ];
L2a.rbc = @(a, c) [a; c];
U2a = L2a \ 0;

% Column 2: [b2; d2] satisfies system (19b) with mu2(L) = 0
L2b = chebop(0, L);
L2b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_R*d - q_R; ...
                      diff(d) - qc_R*b ];
L2b.rbc = @(b, d) [b; d];
U2b = L2b \ 0;

%% 5. Evaluate Jost Solutions at the Interface x = 0
% Reconstruct mu_j(0) = [a_j(0), b_j(0); c_j(0), d_j(0)]
mu1_0 = [U1a{1}(0), U1b{1}(0); ...
         U1a{2}(0), U1b{2}(0)];

mu2_0 = [U2a{1}(0), U2b{1}(0); ...
         U2a{2}(0), U2b{2}(0)];

%% 6. Assemble Scattering Matrix S(lambda)
% S(lambda) = [phi_2(0)]^{-1} * phi_1(0) = [mu_2(0) + I]^{-1} * [mu_1(0) + I]
I = eye(2);
S = (mu2_0 + I) \ (mu1_0 + I);

s11 = S(1,1);
s12 = S(1,2);
s21 = S(2,1);
s22 = S(2,2);

%% 7. Reflection Coefficients (DNNLSE Conventions)
% Defined in Section 2.1, Eq. (8) without artificial signs
rho     = s12 / s22;
rho_hat = s21 / s11;

% Unimodular determinant verification: det(S) = s11*s22 - s12*s21 = 1
det_error = abs(det(S) - 1.0);

%% 8. Display Formatted Results
fprintf('=================================================================\n');
fprintf('  Direct Scattering Transform: Spectral Evaluation Summary\n');
fprintf('=================================================================\n');
fprintf('  Spectral parameter lambda : %+.4f %+.4fi\n', real(lambda), imag(lambda));
fprintf('  Domain truncation length L: %.2f\n', L);
fprintf('-----------------------------------------------------------------\n');
fprintf('  Scattering Matrix Elements:\n');
fprintf('    s11 = %12.8f %+12.8fi\n', real(s11), imag(s11));
fprintf('    s12 = %12.8f %+12.8fi\n', real(s12), imag(s12));
fprintf('    s21 = %12.8f %+12.8fi\n', real(s21), imag(s21));
fprintf('    s22 = %12.8f %+12.8fi\n', real(s22), imag(s22));
fprintf('-----------------------------------------------------------------\n');
fprintf('  Reflection Coefficients:\n');
fprintf('    rho     = %12.8f %+12.8fi\n', real(rho), imag(rho));
fprintf('    rho_hat = %12.8f %+12.8fi\n', real(rho_hat), imag(rho_hat));
fprintf('-----------------------------------------------------------------\n');
fprintf('  Wronskian Identity Check:\n');
fprintf('    |det(S) - 1| = %.4e\n', det_error);
fprintf('=================================================================\n');