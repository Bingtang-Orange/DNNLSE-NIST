clear; clc;
warning('off', 'chebfun:constructor');

L = 20;
lambda = 0.5;
q_full = chebfun(@(x) exp(-x.^2), [-L, L], 'splitting', 'on');

q_L  = chebfun(@(x) q_full(x),        [-L, 0]);
qc_L = chebfun(@(x) conj(q_full(-x)), [-L, 0]);

q_R  = chebfun(@(x) q_full(x),        [0, L]);
qc_R = chebfun(@(x) conj(q_full(-x)), [0, L]);

L1a = chebop(-L, 0);
L1a.op = @(x, a, c) [ diff(a) - q_L*c; ...
                      diff(c) - qc_L*a - 2i*lambda*c - qc_L ];
L1a.lbc = @(a, c) [a; c];
U1a = L1a \ 0;

L1b = chebop(-L, 0);
L1b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_L*d - q_L; ...
                      diff(d) - qc_L*b ];
L1b.lbc = @(b, d) [b; d];
U1b = L1b \ 0;

L2a = chebop(0, L);
L2a.op = @(x, a, c) [ diff(a) - q_R*c; ...
                      diff(c) - qc_R*a - 2i*lambda*c - qc_R ];
L2a.rbc = @(a, c) [a; c];
U2a = L2a \ 0;

L2b = chebop(0, L);
L2b.op = @(x, b, d) [ diff(b) + 2i*lambda*b - q_R*d - q_R; ...
                      diff(d) - qc_R*b ];
L2b.rbc = @(b, d) [b; d];
U2b = L2b \ 0;

mu1_0 = [U1a{1}(0), U1b{1}(0); ...
         U1a{2}(0), U1b{2}(0)];

mu2_0 = [U2a{1}(0), U2b{1}(0); ...
         U2a{2}(0), U2b{2}(0)];

I = eye(2);
S = (mu2_0 + I) \ (mu1_0 + I);

s11 = S(1,1);
s12 = S(1,2);
s21 = S(2,1);
s22 = S(2,2);

rho = s12 / s22;
rho_hat = s21 / s11;
det_error = abs(det(S) - 1.0);

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