% RUN_FIG_RHP2_EVOLUTION
% Simulates multi-time-scale wavepacket dynamics via RHP II solver.
% Reproduces Figures 16-18 (t = 1.0, 100.0, 500.0) with real/imaginary profiles
% and black modulus envelopes +/- |q(x, t)|.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% ---------- 1. Benchmark Potential and Task Configuration ----------
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

% Task specifications: {Time t, Spatial Window [x_min, x_max], Step dx}
tasks = {
    1.0,   [-30, 30],       0.2;
    100.0, [-800, 800],     1.0;
    500.0, [-4000, 4000],   4.0
};

% Solver discretization parameters
L_space    = 20.0;  % Spatial truncation for direct scattering
L_spectral = 8.0;   % Truncation radius of real spectral axis
n_ray      = 24;    % Collocation nodes per diagonal steepest descent ray
n_real     = 64;    % Collocation nodes on real segment Gamma_2

% ---------- 2. Parallel Pool Initialization ----------
pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

num_tasks    = size(tasks, 1);
time_records = zeros(num_tasks, 1);
results      = cell(num_tasks, 1);

% ---------- 3. Main Computation Loop ----------
for k_task = 1:num_tasks
    t_val = tasks{k_task, 1};
    x_lim = tasks{k_task, 2};
    dx    = tasks{k_task, 3};

    x_grid  = (x_lim(1):dx:x_lim(2)).';
    num_pts = length(x_grid);

    fprintf('==================================================================================\n');
    fprintf('>>> Processing Task [%d/%d]: t = %g, Interval = [%g, %g], dx = %g (%d points)\n', ...
            k_task, num_tasks, t_val, x_lim(1), x_lim(2), dx, num_pts);
    fprintf('==================================================================================\n');

    % --- Stage A: Zero-sum condition diagnostic at representative sample points ---
    diag_x = [x_lim(1)*0.75, x_lim(1)*0.25, 0.0, x_lim(2)*0.25, x_lim(2)*0.75];

    fprintf('[Zero-Sum Condition Diagnostic Sampling]\n');
    fprintf('%-10s | %-10s | %-24s | %-14s\n', 'x', 'lambda0', 'q_num(x, t)', 'Zero-Sum Res');
    fprintf('----------------------------------------------------------------------------------\n');
    for j = 1:length(diag_x)
        x_check = diag_x(j);
        [qv_diag, info_diag] = solve_rhp2(x_check, t_val, q0, L_space, L_spectral, n_ray, n_real);
        fprintf('%+10.2f | %+10.4f | %+.6e%+.6ei | %.4e\n', ...
                x_check, info_diag.lambda0, real(qv_diag), imag(qv_diag), info_diag.zero_sum_max);
    end
    fprintf('----------------------------------------------------------------------------------\n');

    % --- Stage B: Parallel dense grid evaluation ---
    q_vals = complex(zeros(num_pts, 1));
    fprintf('Executing parallel grid solve (%d evaluation points) ...\n', num_pts);

    t_start = tic;
    parfor j = 1:num_pts
        [q_vals(j), ~] = solve_rhp2(x_grid(j), t_val, q0, L_space, L_spectral, n_ray, n_real);
    end
    t_elapsed = toc(t_start);
    time_records(k_task) = t_elapsed;

    fprintf('>> Task completed in %.2f s (Average: %.4f s/point).\n\n', ...
            t_elapsed, t_elapsed / num_pts);

    results{k_task}.t    = t_val;
    results{k_task}.x    = x_grid;
    results{k_task}.q    = q_vals;
    results{k_task}.time = t_elapsed;

    q_abs = abs(q_vals);

    % --- Stage C: Multi-Panel Plotting (Figures 16-18 layout) ---
    figure('Color', 'w', 'Position', [100 + 40*k_task, 80 + 30*k_task, 900, 600]);

    % 1. Real part and envelope (Upper panel)
    subplot(2, 1, 1);
    plot(x_grid, real(q_vals), 'r-',  'LineWidth', 1.1, 'DisplayName', '$\mathrm{Re}\,q(x,t)$'); hold on;
    plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.3, 'DisplayName', '$+|q(x,t)|$');
    plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.3, 'DisplayName', '$-|q(x,t)|$'); hold off;
    grid on; box on;
    xlim(x_lim);
    set(gca, 'FontSize', 10, 'LineWidth', 1.0);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\mathrm{Re}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Wavepacket Profile with Envelope ($t = %g$, Time: %.1f s, Points: %d)', ...
          t_val, t_elapsed, num_pts), 'Interpreter', 'latex', 'FontSize', 11);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

    % 2. Imaginary part and envelope (Lower panel)
    subplot(2, 1, 2);
    plot(x_grid, imag(q_vals), 'b-',  'LineWidth', 1.1, 'DisplayName', '$\mathrm{Im}\,q(x,t)$'); hold on;
    plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.3, 'DisplayName', '$+|q(x,t)|$');
    plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.3, 'DisplayName', '$-|q(x,t)|$'); hold off;
    grid on; box on;
    xlim(x_lim);
    set(gca, 'FontSize', 10, 'LineWidth', 1.0);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\mathrm{Im}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Imaginary Component with Envelope ($t = %g$)', t_val), ...
          'Interpreter', 'latex', 'FontSize', 11);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
end

% ---------- 4. Summary Timing Report ----------
fprintf('==================================================================================\n');
fprintf('                           Multi-Time Performance Report                          \n');
fprintf('==================================================================================\n');
fprintf('%-8s | %-16s | %-8s | %-10s | %-12s | %-12s\n', ...
        'Time t', 'Spatial Window', 'Step dx', 'Total Pts', 'Wall Time (s)', 'Time/Pt (s)');
fprintf('----------------------------------------------------------------------------------\n');
for k = 1:num_tasks
    t_v  = tasks{k, 1};
    x_rg = tasks{k, 2};
    d_x  = tasks{k, 3};
    n_p  = length(results{k}.x);
    t_el = time_records(k);
    fprintf('%8g | [%+6g, %+6g] | %8g | %10d | %12.2f | %12.4f\n', ...
            t_v, x_rg(1), x_rg(2), d_x, n_p, t_el, t_el / n_p);
end
fprintf('==================================================================================\n');