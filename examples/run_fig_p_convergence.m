% RUN_FIG_P_CONVERGENCE
% Demonstrates geometric spectral convergence (p-refinement) of RHP IV solver.
% Reproduces Figure 2 in the Response to Reviewers and Section 5.2.3 results.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% 1. Evaluation Target and Potential Setup
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
x_test = -300.0;
t_test = 100.0;

% 2. Adaptive Iteration Parameters
tol_stop = 1e-13;      % Stopping tolerance
n_step   = 4;          % Collocation increment
n_curr   = 8;          % Initial collocation nodes per edge
max_iter = 12;

n_records    = [];
q_records    = [];
diff_records = [];

fprintf('====================================================================\n');
fprintf(' Adaptive Chebyshev Spectral Convergence Test (x = %g, t = %g)\n', x_test, t_test);
fprintf('====================================================================\n');

% Compute baseline solution at n = 8
[q_prev, ~] = solve_rhp4(x_test, t_test, q0, 20.0, 8.0, n_curr, n_curr);
n_records = [n_records; n_curr];
q_records = [q_records; q_prev];

for iter = 1:max_iter
    n_curr = n_curr + n_step;
    [q_curr, info] = solve_rhp4(x_test, t_test, q0, 20.0, 8.0, n_curr, n_curr);
    
    % Successive relative difference (Cauchy convergence metric)
    rel_diff = abs(q_curr - q_prev) / abs(q_curr);
    
    n_records    = [n_records; n_curr];
    q_records    = [q_records; q_curr];
    diff_records = [diff_records; rel_diff];
    
    fprintf('  n = %2d:  q = %+.10e %+.10ei,  Step Diff = %.3e,  Zero-Sum = %.2e\n', ...
            n_curr, real(q_curr), imag(q_curr), rel_diff, info.zero_sum_max);
    
    if rel_diff <= tol_stop
        fprintf('>> Stopping tolerance reached (%.1e). Terminating refinement.\n', tol_stop);
        break;
    end
    
    q_prev = q_curr;
end

% 3. Publication Figure Generation (Figure 2 Layout)
plot_n = n_records(2:end);

fig = figure('Color', 'w', 'Position', [150, 150, 620, 420]);
semilogy(plot_n, diff_records, 'ro-', 'LineWidth', 1.3, 'MarkerSize', 6.5, ...
         'MarkerFaceColor', 'r', 'DisplayName', '$|q_{n} - q_{n-\Delta n}| / |q_{n}|$'); hold on;
yline(eps, 'k:', 'LineWidth', 1.1, 'DisplayName', 'Machine Precision ($\varepsilon_{\mathrm{mach}}$)');
hold off;

grid on; box on;
xlim([min(plot_n) - 2, max(plot_n) + 2]);
ylim([1e-16, 1e-1]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('Chebyshev Collocation Nodes per Contour ($n$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Successive Relative Difference', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Geometric Spectral Convergence of RHP IV Solver ($t = %g$)', t_test), ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);