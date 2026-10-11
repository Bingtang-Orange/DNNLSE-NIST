% Comparison of numerical solution and asymptotic envelope at t = 30.

clear; clc; close all;

%% Load data and setup
data_file = 't30_rhp4_data.mat';
if ~exist(data_file, 'file')
    error('Data file %s not found. Please ensure it exists in the working path.', data_file);
end
load(data_file);

q_abs = abs(q_vals);

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

%% Compute asymptotic envelope
x_asymp = (-250:10:250).';
fprintf('Computing asymptotic envelope (%d points)...\n', length(x_asymp));
[~, mod_asymp] = eval_dnnlse_asymptotics(x_asymp, t_val, q0);
fprintf('Done.\n');

%% Plot solution and envelope
figure('Color', 'w', 'Position', [100, 100, 1050, 700]);

subplot(2, 1, 1);
plot(x_grid, real(q_vals), 'r-', 'LineWidth', 0.5, 'DisplayName', 'Re q'); hold on;
plot(x_grid,  q_abs, 'k-',  'LineWidth', 1.0, 'DisplayName', 'Numerical |q|');
plot(x_grid, -q_abs, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(x_asymp,  mod_asymp, 'ko', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic |q|');
plot(x_asymp, -mod_asymp, 'ko', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
     'MarkerFaceColor', 'none', 'HandleVisibility', 'off');
grid on; box on;
xlim([-250, 250]);
ylim([-0.16, 0.16]);
ylabel('Re q(x, 30)', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);

subplot(2, 1, 2);
plot(x_grid, imag(q_vals), 'b-', 'LineWidth', 0.5, 'DisplayName', 'Im q'); hold on;
plot(x_grid,  q_abs, 'k-',  'LineWidth', 1.0, 'DisplayName', 'Numerical |q|');
plot(x_grid, -q_abs, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(x_asymp,  mod_asymp, 'ko', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic |q|');
plot(x_asymp, -mod_asymp, 'ko', 'MarkerSize', 4.5, 'LineWidth', 1.0, ...
     'MarkerFaceColor', 'none', 'HandleVisibility', 'off');
grid on; box on;
xlim([-250, 250]);
ylim([-0.16, 0.16]);
xlabel('x', 'FontSize', 12);
ylabel('Im q(x, 30)', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);