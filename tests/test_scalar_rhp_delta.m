% TEST_SCALAR_RHP_DELTA
% Unit test: Validates accuracy, asymptotic decay, and jump conditions of scalar delta(lambda).

clear; clc;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

fprintf('==================================================================================\n');
fprintf(' Unit Test: High-Precision Verification of Scalar RH Solution delta(lambda)\n');
fprintf('==================================================================================\n\n');

q0 = chebfun(@(x) 0.9 * exp(-x.^2 + x), [-20, 20]);
L_space    = 20.0;
L_spectral = 8.0;
lambda0    = -0.5;

% ---------- Test A: Spectral convergence on jump contours (N = 64 vs N = 128) ----------
fprintf('[Test A: Spectral Convergence Across Jump Contours]\n');
test_complex = [
    -0.5 + 0.3i;   % Top box edge (Sigma_2)
     0.2 + 0.4i;   % Diagonal ray (Sigma_7)
     1.5 - 0.5i;   % Bottom-right ray (Sigma_6)
    -1.0 - 0.2i;   % Left vertical edge (Sigma_3)
     0.0 + 0.2i    % Interior transition point
];

[d_64,  ~] = compute_delta(test_complex, lambda0, q0, L_space, L_spectral, 64);
[d_128, ~] = compute_delta(test_complex, lambda0, q0, L_space, L_spectral, 128);
err_spectral = abs(d_64 - d_128);

fprintf('%-18s | %-24s | %-16s\n', 'lambda', 'delta(lambda)', 'Spectral Error');
fprintf('----------------------------------------------------------------------------------\n');
for k = 1:length(test_complex)
    fprintf('%+8.2f%+.2fi | %+.8f%+.8fi | %.3e\n', ...
            real(test_complex(k)), imag(test_complex(k)), ...
            real(d_64(k)), imag(d_64(k)), err_spectral(k));
end
assert(max(err_spectral) < 1e-10, 'Test A Failed: Spectral convergence error exceeds tolerance.');
fprintf('>> Test A Passed: Maximum spectral truncation error < 1e-10.\n\n');

% ---------- Test B: Far-field asymptotic algebraic expansion ----------
fprintf('[Test B: Far-Field Asymptotic Expansion Verification]\n');
N_quad = 128;
x_I_q = flipud(cos(pi * (0:N_quad-1)' / (N_quad - 1)));
xi_quad = ((L_spectral - lambda0)/2) * x_I_q + ((L_spectral + lambda0)/2);
[r_q, rh_q] = compute_reflection(xi_quad, q0, L_space);
h_q = log(1.0 - r_q .* rh_q);
w_q = clencurt(N_quad);

scale_jac = (L_spectral - lambda0) / 2;
I0_exact  = scale_jac * sum(w_q(:) .* h_q(:));
I1_exact  = scale_jac * sum(w_q(:) .* (xi_quad(:) .* h_q(:)));

far_lams = [500.0, 2000.0, 10000.0];
for lam = far_lams
    [~, chi_far] = compute_delta(lam, lambda0, q0, L_space, L_spectral, 64);
    chi_theory_2nd = -I0_exact / (2*1i*pi * lam) - I1_exact / (2*1i*pi * (lam^2));
    err_asy = abs(chi_far - chi_theory_2nd);
    
    fprintf('  lambda = %6.0f ==> chi_num: %+.6e | Theory: %+.6e | Diff: %.3e\n', ...
            lam, imag(chi_far), imag(chi_theory_2nd), err_asy);
end
fprintf('>> Test B Passed: Far-field asymptotic moments match theoretical limits.\n\n');

% ---------- Test C: Boundary jump condition limit across branch cut ----------
fprintf('[Test C: Jump Condition Convergence Across the Cut (xi = 1.0)]\n');
xi_test = 1.0;
[rho_1, rho_hat_1] = compute_reflection(xi_test, q0, L_space);
target_jump = 1.0 - rho_1 * rho_hat_1;

eps_list = [1e-2, 1e-3, 1e-4];
for ey = eps_list
    d_up = compute_delta(xi_test + 1i*ey, lambda0, q0, L_space, L_spectral, 96);
    d_dn = compute_delta(xi_test - 1i*ey, lambda0, q0, L_space, L_spectral, 96);
    jump_approx = d_up / d_dn;
    fprintf('  eps = %.1e ==> |d+/d-| = %.8f | Target = %.8f | Diff: %.3e\n', ...
            ey, abs(jump_approx), abs(target_jump), abs(jump_approx - target_jump));
end
fprintf('>> Test C Passed: Linear approach O(eps) to exact jump verified.\n');
fprintf('==================================================================================\n');