clear; clc; close all;

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
x_test = -300.0;
t_test = 100.0;

L_space    = 8.0;
L_spectral = 8.0;

n_ref  = 120;             
n_list = (8:6:68).';       
num_n  = length(n_list);

fprintf('====================================================================\n');
fprintf('  Chebyshev Spectral Convergence Test: (x = %g, t = %g)\n', x_test, t_test);
fprintf('====================================================================\n');

fprintf('Computing reference benchmark with n_ref = %d...\n', n_ref);
t_ref_start = tic;
[q_ref, info_ref] = solve_rhp4(x_test, t_test, q0, L_space, L_spectral, n_ref, n_ref);
t_ref_cost = toc(t_ref_start);

fprintf('  q_ref = %+.15e %+.15ei (computed in %.2f s)\n', ...
        real(q_ref), imag(q_ref), t_ref_cost);
fprintf('  Reference Zero-Sum Residual: %.2e\n\n', info_ref.zero_sum_max);

q_records   = complex(zeros(num_n, 1));
err_records = zeros(num_n, 1);
zs_records  = zeros(num_n, 1);

fprintf('%-6s | %-32s | %-14s | %-12s\n', 'n', 'q_n(x, t)', 'Rel Error \epsilon_n', 'Zero-Sum Res');
fprintf('--------------------------------------------------------------------\n');

for k = 1:num_n
    n_curr = n_list(k);
    [q_curr, info] = solve_rhp4(x_test, t_test, q0, L_space, L_spectral, n_curr, n_curr);
    
    rel_err = abs(q_curr - q_ref) / abs(q_ref);
    
    q_records(k)   = q_curr;
    err_records(k) = rel_err;
    zs_records(k)  = info.zero_sum_max;
    
    fprintf('  %2d   | %+.10e %+.10ei |   %.4e   |   %.2e\n', ...
            n_curr, real(q_curr), imag(q_curr), rel_err, info.zero_sum_max);
end
fprintf('--------------------------------------------------------------------\n\n');

figure('Color', 'w', 'Position', [150, 150, 650, 450]);

semilogy(n_list, err_records, 'ro-', 'LineWidth', 1.3, 'MarkerSize', 7, ...
         'MarkerFaceColor', 'r', 'DisplayName', '$|q_{n} - q_{\mathrm{ref}}| / |q_{\mathrm{ref}}|$');

grid on; box on;
xlim([min(n_list) - 2, max(n_list) + 2]);
ylim([1e-12, 1e-2]); 
xlabel('Chebyshev Collocation Points per Contour ($n$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Relative Error $\varepsilon_n$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Spectral Convergence of RHP IV Solver ($t = %g$)', t_test), 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);