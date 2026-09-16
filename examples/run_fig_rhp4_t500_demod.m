% RUN_FIG_RHP4_T500_DEMOD
% Reconstructs the slow-varying demodulated wavepacket for DNNLSE at t = 500.
% Reproduces Figure 18 by removing rapid oscillations via carrier demodulation
% factor exp(+i*x^2/(4t)).

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% 1. Configuration and Parameters
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);
t_val = 500.0;

x_min   = -5000.0;
x_max   =  5000.0;
dx      = 50.0; 
x_grid  = (x_min:dx:x_max).';
num_pts = length(x_grid);

fprintf('==================================================================================\n');
fprintf(' RHP IV Large-Time Evaluation and Carrier Demodulation (t = %g)\n', t_val);
fprintf(' Spatial range: [%d, %d], Step size dx = %g (%d evaluation points)\n', ...
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

% 3. Parallel Solve Across Spatial Grid
t_start = tic;
parfor j = 1:num_pts
    [qv, info] = solve_rhp4(x_grid(j), t_val, q0, 20.0, 8.0, 24, 20);
    q_vals(j)       = qv;
    zerosum_vals(j) = info.zero_sum_max;
    send(data_queue, j);
end
t_cost = toc(t_start);

fprintf('>> Solve completed in %.2f s (Average: %.4f s/point).\n', t_cost, t_cost / num_pts);
fprintf('>> Maximum zero-sum residual: %.3e\n\n', max(zerosum_vals));

% 4. Demodulate Carrier Phase: q_demod = q * exp(+i * x^2 / (4*t))
fast_phase = exp(+1i * (x_grid.^2) / (4 * t_val));
q_demod    = q_vals .* fast_phase;
q_abs      = abs(q_vals);

% 5. Demodulated Waveform Visualization (Figure 18 Layout)
fig = figure('Color', 'w', 'Position', [150, 150, 850, 480]);

plot(x_grid, real(q_demod), 'r-',  'LineWidth', 1.2, 'DisplayName', '$\mathrm{Re}[q\,\mathrm{e}^{+\mathrm{i}x^2/(4t)}]$'); hold on;
plot(x_grid, imag(q_demod), 'b--', 'LineWidth', 1.2, 'DisplayName', '$\mathrm{Im}[q\,\mathrm{e}^{+\mathrm{i}x^2/(4t)}]$');
plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.3, 'DisplayName', 'Envelope $|q|$');
plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.3, 'HandleVisibility', 'off'); hold off;

grid on; box on;
xlim([x_min, x_max]);
set(gca, 'FontSize', 11, 'LineWidth', 1.0);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Demodulated Amplitude', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('Demodulated Waveform and Amplitude Envelope ($t = %g$)', t_val), ...
      'Interpreter', 'latex', 'FontSize', 13);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);

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