% Spectral convergence (p-refinement) of RHP IV solver with respect to collocation points n.

clear; clc; close all;

%% Parameters and setup
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
x_test = -300.0;
t_test = 100.0;

tol_stop = 1e-13;
n_step   = 4;
n_curr   = 8;
max_iter = 12;

n_records = [];
q_records = [];
diff_records = [];

fprintf('====================================================================\n');
fprintf('  Chebyshev Spectral Convergence Test: (x = %g, t = %g)\n', x_test, t_test);
fprintf('====================================================================\n');

[q_prev, ~] = solve_rhp4(x_test, t_test, q0, 8.0, 8.0, n_curr, n_curr);
n_records = [n_records; n_curr];
q_records = [q_records; q_prev];

%% Convergence iteration
for iter = 1:max_iter
    n_curr = n_curr + n_step;
    [q_curr, info] = solve_rhp4(x_test, t_test, q0, 8.0, 8.0, n_curr, n_curr);
    
    rel_diff = abs(q_curr - q_prev) / abs(q_curr);
    
    n_records = [n_records; n_curr];
    q_records = [q_records; q_curr];
    diff_records = [diff_records; rel_diff];
    
    fprintf('  n = %2d:  q = %+.10e %+.10ei,  Rel Diff = %.3e,  Zero-Sum = %.2e\n', ...
            n_curr, real(q_curr), imag(q_curr), rel_diff, info.zero_sum_max);
    
    if rel_diff <= tol_stop
        fprintf('Tolerance %.1e reached. Stopping iteration.\n', tol_stop);
        break;
    end
    
    q_prev = q_curr;
end

%% Plot convergence
plot_n = n_records(2:end);

figure('Color', 'w', 'Position', [150, 150, 650, 450]);
semilogy(plot_n, diff_records, 'ro-', 'LineWidth', 1.3, 'MarkerSize', 7, ...
         'MarkerFaceColor', 'r', 'DisplayName', '$|q_{n} - q_{n-\Delta n}| / |q_{n}|$'); hold on;
yline(eps, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Machine Precision ($\varepsilon$)');

grid on; box on;
xlim([min(plot_n)-2, max(plot_n)+2]);
ylim([1e-16, 1e-1]);
xlabel('Chebyshev Collocation Points per Contour ($n$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Successive Relative Difference', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Spectral Convergence of RHP IV Solver ($t = %g$)', t_test), 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);