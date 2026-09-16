% RUN_FIG_RHP4_T100_ENVELOPE
% Simulates long-time wavepacket dynamics via RHP IV at t = 100.
% Computes real/imaginary wavepacket profiles, modulus envelope +/- |q(x, t)|,
% and spatial distribution of the zero-sum compatibility residual.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% 1. Configuration and Parameters
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
t_val = 100.0;

x_min   = -1000.0;
x_max   =  900.0;
dx      = 2.0;
x_grid  = (x_min:dx:x_max).';
num_pts = length(x_grid);

fprintf('==================================================================================\n');
fprintf(' RHP IV Wavepacket Dynamics and Zero-Sum Diagnostic (t = %g)\n', t_val);
fprintf(' Spatial domain: [%d, %d], Step size dx = %g (%d evaluation points)\n', ...
        x_min, x_max, dx, num_pts);
fprintf('==================================================================================\n');

% 2. Parallel Pool Initialization
pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

q_vals       = complex(zeros(num_pts, 1));
zerosum_vals = zeros(num_pts, 1);
data_queue   = parallel.pool.DataQueue;
afterEach(data_queue, @(idx) progress_counter(num_pts));

% 3. Parallel Dense Grid Computation
t_start = tic;
parfor j = 1:num_pts
    [qv, info] = solve_rhp4(x_grid(j), t_val, q0, 20.0, 8.0, 24, 20);
    q_vals(j)       = qv;
    zerosum_vals(j) = info.zero_sum_max;
    send(data_queue, j);
end
t_cost = toc(t_start);

fprintf('>> Solve completed in %.2f s (Average: %.4f s/point).\n', t_cost, t_cost / num_pts);
fprintf('>> Maximum zero-sum residual across entire domain: %.3e\n\n', max(zerosum_vals));

q_abs = abs(q_vals);

% 4. Three-Panel Diagnostic Visualization
figure('Color', 'w', 'Position', [100, 80, 950, 750]);

% Panel (a): Real component and envelope
subplot(3, 1, 1);
plot(x_grid, real(q_vals), 'r-',  'LineWidth', 1.0, 'DisplayName', '$\mathrm{Re}\,q(x, t)$'); hold on;
plot(x_grid,  q_abs,       'k--', 'LineWidth', 1.1, 'DisplayName', 'Envelope $\pm|q(x, t)|$');
plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.1, 'HandleVisibility', 'off'); hold off;
grid on; box on;
xlim([x_min, x_max]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\mathrm{Re}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Real Component and Envelopes at $t = %g$ (Wall time: %.1f s)', t_val, t_cost), ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

% Panel (b): Imaginary component and envelope
subplot(3, 1, 2);
plot(x_grid, imag(q_vals), 'b--', 'LineWidth', 1.0, 'DisplayName', '$\mathrm{Im}\,q(x, t)$'); hold on;
plot(x_grid,  q_abs,       'k--', 'LineWidth', 1.1, 'DisplayName', 'Envelope $\pm|q(x, t)|$');
plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.1, 'HandleVisibility', 'off'); hold off;
grid on; box on;
xlim([x_min, x_max]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\mathrm{Im}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title('Imaginary Component and Envelopes', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

% Panel (c): Maximum zero-sum junction residual across the domain
subplot(3, 1, 3);
semilogy(x_grid, zerosum_vals, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 1.2, ...
         'DisplayName', 'Max Zero-Sum Residual'); hold on;
yline(1e-10, 'r:', 'LineWidth', 1.1, 'DisplayName', 'Threshold ($10^{-10}$)'); hold off;
grid on; box on;
xlim([x_min, x_max]);
ylim([1e-15, 1e-4]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Zero-Sum Residual', 'Interpreter', 'latex', 'FontSize', 11);
title('Five-Node Zero-Sum Compatibility Metric (RHP IV Health Check)', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

% Progress Counter Callback
function progress_counter(total_pts)
    persistent current_count
    if isempty(current_count)
        current_count = 0;
    end
    current_count = current_count + 1;
    fprintf('\r[Progress] Completed: %4d / %4d (%.1f%%)', ...
            current_count, total_pts, 100 * (current_count / total_pts));
    if current_count == total_pts
        current_count = [];
        fprintf('\n');
    end
end