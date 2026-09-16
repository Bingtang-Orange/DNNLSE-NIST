% RUN_FIG_CROSS_VALIDATION
% Evaluates cross-stage consistency between RHP II and RHP IV at t = 1.0.
% Computes 151-point smooth profiles and pointwise absolute discrepancies.
% Reproduces Figure 1 in the Response to Reviewers and Section 5.2.3 results.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% 1. Benchmark Potential and High-Density Grid
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
t_val = 1.0;

x_min   = -15.0;
x_max   =  15.0;
dx      =  0.2; 
x_grid  = (x_min:dx:x_max).';
num_pts = length(x_grid);

fprintf('====================================================================\n');
fprintf(' Cross-Validation: RHP II vs RHP IV (t = 1.0, Grid Points: %d)\n', num_pts);
fprintf('====================================================================\n');

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

q_rhp2 = complex(zeros(num_pts, 1));
q_rhp4 = complex(zeros(num_pts, 1));

t_start = tic;
parfor j = 1:num_pts
    xev = x_grid(j);
    
    % RHP II solver: n_ray = 32, n_real = 64
    [q2, ~] = solve_rhp2(xev, t_val, q0, 20.0, 8.0, 32, 64);
    q_rhp2(j) = q2;
    
    % RHP IV solver: explicit geometric configuration r_box = 0.25, d_ray = 2.50
    [q4, ~] = solve_rhp4(xev, t_val, q0, 20.0, 8.0, 24, 32, 0.25, 2.50);
    q_rhp4(j) = q4;
end
t_cost = toc(t_start);

abs_diff = abs(q_rhp2 - q_rhp4);
fprintf('>> Solve completed in %.2f s (Average: %.4f s/point).\n', t_cost, t_cost / num_pts);
fprintf('>> Maximum pointwise absolute discrepancy: %.3e\n\n', max(abs_diff));

% 2. Two-Panel Diagnostic Figure (Figure 1 Layout)
fig = figure('Color', 'w', 'Position', [100, 150, 950, 400]);

% Panel (a): Reconstructed Potential Profiles
subplot(1, 2, 1);
plot(x_grid, real(q_rhp2), 'r-',  'LineWidth', 1.3, 'DisplayName', '$\mathrm{Re}\,q$ (RHP II)'); hold on;
plot(x_grid, imag(q_rhp2), 'b--', 'LineWidth', 1.3, 'DisplayName', '$\mathrm{Im}\,q$ (RHP II)');

% Sparse markers for RHP IV to ensure visual clarity
idx_mark = 1:5:num_pts;
plot(x_grid(idx_mark), real(q_rhp4(idx_mark)), 'ko', 'MarkerSize', 4.5, ...
     'LineWidth', 1.0, 'DisplayName', '$\mathrm{Re}\,q$ (RHP IV)');
plot(x_grid(idx_mark), imag(q_rhp4(idx_mark)), 'ks', 'MarkerSize', 4.5, ...
     'LineWidth', 1.0, 'DisplayName', '$\mathrm{Im}\,q$ (RHP IV)');
hold off;

grid on; box on;
xlim([x_min, x_max]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$q(x, 1)$', 'Interpreter', 'latex', 'FontSize', 12);
title('(a) Solution Profiles at $t = 1$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);

% Panel (b): Absolute Pointwise Discrepancy (Semilog)
subplot(1, 2, 2);
semilogy(x_grid, abs_diff, 'm-', 'LineWidth', 1.3, ...
         'DisplayName', '$|q^{\mathrm{RHP\,II}} - q^{\mathrm{RHP\,IV}}|$'); hold on;
semilogy(x_grid(idx_mark), abs_diff(idx_mark), 'm.', 'MarkerSize', 9, 'HandleVisibility', 'off');
hold off;

grid on; box on;
xlim([x_min, x_max]);
ylim([1e-9, 1e-5]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Absolute Difference', 'FontSize', 11);
title('(b) Pointwise Discrepancy', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);