clear; clc; close all;

mat_file = 't30_rhp4_data.mat';
if ~exist(mat_file, 'file')
    error('Data file %s not found. Please ensure it is in the current working directory.', mat_file);
end
load(mat_file);

idx_core = (x_grid >= -50.0) & (x_grid <= 50.0);
x_core   = x_grid(idx_core);
q_core   = q_vals(idx_core);

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

dx_asymp  = 1.0; 
x_asymp   = (-50.0:dx_asymp:50.0).';
[q_asymp, ~] = eval_dnnlse_asymptotics(x_asymp, t_val, q0);

fig = figure('Name', 'Core Region Zoom', 'Color', 'w', ...
             'Position', [150, 150, 1050, 650]);

subplot(2, 1, 1);
plot(x_core, real(q_core), 'r-', 'LineWidth', 1.1, 'DisplayName', 'Numerical Re q'); hold on;
plot(x_asymp, real(q_asymp), 'ko', 'MarkerSize', 4.0, 'LineWidth', 0.8, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic Re q');
grid on; box on;
xlim([-50, 50]);
ylim([-0.15, 0.15]);
ylabel('$\mathrm{Re}\,q(x, 30)$', 'Interpreter', 'latex', 'FontSize', 12);
title('Core Region Zoom: $x \in [-50, 50]$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');

subplot(2, 1, 2);
plot(x_core, imag(q_core), 'b-', 'LineWidth', 1.1, 'DisplayName', 'Numerical Im q'); hold on;
plot(x_asymp, imag(q_asymp), 'ko', 'MarkerSize', 4.0, 'LineWidth', 0.8, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic Im q');
grid on; box on;
xlim([-50, 50]);
ylim([-0.15, 0.15]);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\mathrm{Im}\,q(x, 30)$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');