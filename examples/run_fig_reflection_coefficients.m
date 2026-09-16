% RUN_FIG_REFLECTION_COEFFICIENTS
% Computes and plots the real part, imaginary part, and modulus of the reflection
% coefficients rho(lambda) and rho_hat(lambda) for Gaussian initial data.
% Reproduces Figure for Section 3.3.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

%% 1. Simulation Parameters
L = 40.0;                        % Truncation length ensuring machine precision
q0 = @(x) exp(-x.^2);            % Gaussian initial potential

% Spectral grid on the real axis
lambda_vec = -2.0:0.05:2.0;
num_lambda = length(lambda_vec);

fprintf('Computing reflection coefficients across lambda in [%.1f, %.1f]...\n', ...
        lambda_vec(1), lambda_vec(end));

%% 2. Vectorized Direct Scattering Evaluation
[rho_arr, rho_hat_arr] = compute_reflection(lambda_vec, q0, L);

%% 3. Generate Paper Figure (Two-Panel Layout)
fig = figure('Color', 'w', 'Position', [150, 150, 750, 600]);

% Panel (a): Left-incident reflection coefficient rho(lambda)
subplot(2, 1, 1);
plot(lambda_vec, real(rho_arr), 'b-', 'LineWidth', 1.5); hold on;
plot(lambda_vec, imag(rho_arr), 'r--', 'LineWidth', 1.5);
plot(lambda_vec, abs(rho_arr),  'k-', 'LineWidth', 2.0);
hold off;
grid on; box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.1);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$\rho(\lambda)$', 'Interpreter', 'latex', 'FontSize', 13);
title('Reflection Coefficient $\rho(\lambda)$ (Gaussian Potential)', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend({'$\mathrm{Re}\,\rho$', '$\mathrm{Im}\,\rho$', '$|\rho|$'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
xlim([-2.0, 2.0]);

% Panel (b): Nonlocal reflection coefficient rho_hat(lambda)
subplot(2, 1, 2);
plot(lambda_vec, real(rho_hat_arr), 'b-', 'LineWidth', 1.5); hold on;
plot(lambda_vec, imag(rho_hat_arr), 'r--', 'LineWidth', 1.5);
plot(lambda_vec, abs(rho_hat_arr),  'k-', 'LineWidth', 2.0);
hold off;
grid on; box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.1);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$\hat{\rho}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 13);
title('Nonlocal Reflection Coefficient $\hat{\rho}(\lambda)$', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend({'$\mathrm{Re}\,\hat{\rho}$', '$\mathrm{Im}\,\hat{\rho}$', '$|\hat{\rho}|$'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
xlim([-2.0, 2.0]);

%% 4. Command-Window Diagnostics Summary
fprintf('\n================== Summary of Spectral Data ==================\n');
fprintf('Grid size      : %d points (step size = %.2f)\n', num_lambda, lambda_vec(2)-lambda_vec(1));
fprintf('Peak |rho|     : %.6e at lambda = %.2f\n', max(abs(rho_arr)), lambda_vec(abs(rho_arr) == max(abs(rho_arr))));
fprintf('Peak |rho_hat| : %.6e at lambda = %.2f\n', max(abs(rho_hat_arr)), lambda_vec(abs(rho_hat_arr) == max(abs(rho_hat_arr))));
fprintf('==============================================================\n');