% Numerical inverse scattering reconstruction using RHP IV at t = 30.

clear; clc; close all;

%% Parameters and spatial grid
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

t_val = 30.0;
x_min = -250.0;
x_max =  250.0;

x_extreme     = max(abs([x_min, x_max]));
k_max         = x_extreme / (2 * t_val);
lambda_min    = 2 * pi / k_max;
pts_per_cycle = 6.0;
dx_theory     = lambda_min / pts_per_cycle;
dx            = round(dx_theory, 2);

x_grid  = (x_min:dx:x_max).';
num_pts = length(x_grid);

fprintf('==================================================================================\n');
fprintf('  DNNLSE RHP IV Solver (t = %g)\n', t_val);
fprintf('  Domain: [%.1f, %.1f], dx = %.2f (%d points)\n', x_min, x_max, dx, num_pts);
fprintf('==================================================================================\n');

%% Start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('Processes');
end
num_workers = pool.NumWorkers;

%% Setup progress queue
q_vals        = complex(zeros(num_pts, 1));
zerosum_vals  = zeros(num_pts, 1);
data_queue    = parallel.pool.DataQueue;
afterEach(data_queue, @(idx) progress_counter(num_pts));

%% Parallel RHP IV solve
t_start = tic;
parfor j = 1:num_pts
    [qv, info] = solve_rhp4(x_grid(j), t_val, q0, 8.0, 8.0, 24, 20);
    q_vals(j)       = qv;
    zerosum_vals(j) = info.zero_sum_max;
    send(data_queue, j);
end
t_cost = toc(t_start);

fprintf('Computation finished in %.2f s (%.4f s/pt).\n', t_cost, t_cost / num_pts);
fprintf('Maximum zero-sum residual: %.3e\n\n', max(zerosum_vals));

%% Plot solution and envelope
q_abs = abs(q_vals);

figure('Color', 'w', 'Position', [100, 100, 1050, 700]);

subplot(2, 1, 1);
plot(x_grid, real(q_vals), 'r-',  'LineWidth', 0.6, 'DisplayName', 'Re q'); hold on;
plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.0, 'DisplayName', '|q|');
plot(x_grid, -q_abs,       'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
grid on; box on;
xlim([x_min, x_max]);
ylabel('Re q(x, 30)', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);

subplot(2, 1, 2);
plot(x_grid, imag(q_vals), 'b--', 'LineWidth', 0.8, 'DisplayName', 'Im q'); hold on;
plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.0, 'DisplayName', '|q|');
plot(x_grid, -q_abs,       'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
grid on; box on;
xlim([x_min, x_max]);
xlabel('x', 'FontSize', 12);
ylabel('Im q(x, 30)', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8);

%% Helper functions
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