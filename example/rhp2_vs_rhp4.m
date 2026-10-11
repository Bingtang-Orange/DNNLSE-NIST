% Cross-validation between independent solvers for RHP II and RHP IV at t = 1.

clear; clc; close all;

%% Parameters and spatial grid
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
t_val = 1.0;

x_min = -15.0;
x_max =  15.0;
dx    =  0.2; 
x_grid = (x_min:dx:x_max).';
num_pts = length(x_grid);

%% Parallel solve for RHP II and RHP IV
pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

q_rhp2 = complex(zeros(num_pts, 1));
q_rhp4 = complex(zeros(num_pts, 1));

fprintf('Computing cross-validation for t = 1.0 (%d points)...\n', num_pts);
t_start = tic;
parfor j = 1:num_pts
    xev = x_grid(j);
    
    [q2, ~] = solve_rhp2(xev, t_val, q0, 8.0, 8.0, 32, 64);
    q_rhp2(j) = q2;
    
    [q4, ~] = solve_rhp4(xev, t_val, q0, 8.0, 8.0, 24, 32, 0.25, 2.50);
    q_rhp4(j) = q4;
end
t_cost = toc(t_start);

abs_diff = abs(q_rhp2 - q_rhp4);
fprintf('Done in %.2f s. Maximum absolute discrepancy: %.3e\n\n', t_cost, max(abs_diff));

%% Plot comparison and discrepancy
fig = figure('Color', 'w', 'Position', [100, 150, 1050, 430]);

subplot(1, 2, 1);
plot(x_grid, real(q_rhp2), 'r-',  'LineWidth', 1.3, 'DisplayName', 'Re q (RHP II)'); hold on;
plot(x_grid, imag(q_rhp2), 'b--', 'LineWidth', 1.3, 'DisplayName', 'Im q (RHP II)');

idx_mark = 1:5:num_pts;
plot(x_grid(idx_mark), real(q_rhp4(idx_mark)), 'ko', 'MarkerSize', 4.5, ...
     'LineWidth', 1.0, 'DisplayName', 'Re q (RHP IV)');
plot(x_grid(idx_mark), imag(q_rhp4(idx_mark)), 'ks', 'MarkerSize', 4.5, ...
     'LineWidth', 1.0, 'DisplayName', 'Im q (RHP IV)');

grid on; box on;
xlim([x_min, x_max]);
xlabel('x', 'FontSize', 11);
ylabel('q(x, 1)', 'FontSize', 11);
title('(a) Solution Profiles at t = 1', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9);

subplot(1, 2, 2);
semilogy(x_grid, abs_diff, 'm-', 'LineWidth', 1.3, ...
         'DisplayName', '|q^{\mathrm{RHP\,II}} - q^{\mathrm{RHP\,IV}}|'); hold on;
semilogy(x_grid(idx_mark), abs_diff(idx_mark), 'm.', 'MarkerSize', 10, 'HandleVisibility', 'off');

grid on; box on;
xlim([x_min, x_max]);
ylim([1e-9, 1e-5]);
xlabel('x', 'FontSize', 11);
ylabel('Absolute Difference', 'FontSize', 11);
title('(b) Pointwise Discrepancy', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);