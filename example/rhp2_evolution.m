clear; clc; close all;

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

tasks = {
    1.0,   [-30, 30],       0.2;
    100.0, [-800, 800],     1.0;
    500.0, [-4000, 4000],   4.0
};

L_space    = 8.0;
L_spectral = 8.0;
n_ray      = 24;
n_real     = 64;

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

time_records = zeros(size(tasks, 1), 1);
results = cell(size(tasks, 1), 1);

for k_task = 1:size(tasks, 1)
    t_val = tasks{k_task, 1};
    x_lim = tasks{k_task, 2};
    dx    = tasks{k_task, 3};
    
    x_grid  = (x_lim(1):dx:x_lim(2)).';
    num_pts = length(x_grid);
    
    fprintf('==================================================================================\n');
    fprintf('Task [%d/%d]: t = %g, domain = [%g, %g], dx = %g (%d points)\n', ...
            k_task, size(tasks, 1), t_val, x_lim(1), x_lim(2), dx, num_pts);
    fprintf('==================================================================================\n');
    
    diag_x = [x_lim(1)*0.75, x_lim(1)*0.25, 0.0, x_lim(2)*0.25, x_lim(2)*0.75];
    
    fprintf('Zero-sum residual diagnostic:\n');
    fprintf('%-10s | %-10s | %-24s | %-12s\n', 'x', 'lambda0', 'q(x, t)', 'Zero-Sum Res');
    fprintf('----------------------------------------------------------------------------------\n');
    for j = 1:length(diag_x)
        x_check = diag_x(j);
        [qv_diag, info_diag] = solve_rhp2(x_check, t_val, q0, L_space, L_spectral, n_ray, n_real);
        fprintf('%+10.2f | %+10.4f | %+.6e%+.6ei | %.3e\n', ...
                x_check, info_diag.lambda0, real(qv_diag), imag(qv_diag), info_diag.zero_sum_max);
    end
    fprintf('----------------------------------------------------------------------------------\n');
    
    q_vals = complex(zeros(num_pts, 1));
    fprintf('Solving in parallel (%d points)...\n', num_pts);
    
    t_start = tic;
    parfor j = 1:num_pts
        [q_vals(j), ~] = solve_rhp2(x_grid(j), t_val, q0, L_space, L_spectral, n_ray, n_real);
    end
    t_elapsed = toc(t_start);
    time_records(k_task) = t_elapsed;
    
    fprintf('Done in %.2f s (%.4f s/point).\n\n', t_elapsed, t_elapsed / num_pts);
    
    results{k_task}.t = t_val;
    results{k_task}.x = x_grid;
    results{k_task}.q = q_vals;
    results{k_task}.time = t_elapsed;
    
    q_abs = abs(q_vals);
    
    figure('Color', 'w', 'Position', [100 + 40*k_task, 80 + 30*k_task, 980, 640]);
    
    subplot(2, 1, 1);
    plot(x_grid, real(q_vals), 'r-',  'LineWidth', 1.0, 'DisplayName', '$\mathrm{Re}\,q(x,t)$'); hold on;
    plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.2, 'DisplayName', '$|q(x,t)|$');
    plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.2, 'DisplayName', '$-|q(x,t)|$');
    grid on; box on;
    xlim(x_lim);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$\mathrm{Re}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('$t = %g$', t_val), 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    
    subplot(2, 1, 2);
    plot(x_grid, imag(q_vals), 'g-',  'LineWidth', 1.0, 'DisplayName', '$\mathrm{Im}\,q(x,t)$'); hold on;
    plot(x_grid,  q_abs,       'k-',  'LineWidth', 1.2, 'DisplayName', '$|q(x,t)|$');
    plot(x_grid, -q_abs,       'k--', 'LineWidth', 1.2, 'DisplayName', '$-|q(x,t)|$');
    grid on; box on;
    xlim(x_lim);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$\mathrm{Im}\,q(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('$t = %g$', t_val), 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
end

fprintf('==================================================================================\n');
fprintf('                           Benchmark Summary                                      \n');
fprintf('==================================================================================\n');
fprintf('%-8s | %-16s | %-8s | %-10s | %-12s | %-12s\n', ...
        't', 'Domain', 'dx', 'Points', 'Time (s)', 'Time/pt (s)');
fprintf('----------------------------------------------------------------------------------\n');
for k = 1:size(tasks, 1)
    t_v  = tasks{k, 1};
    x_rg = tasks{k, 2};
    d_x  = tasks{k, 3};
    n_p  = length(results{k}.x);
    t_el = time_records(k);
    fprintf('%8g | [%+6g, %+6g] | %8g | %10d | %12.2f | %12.4f\n', ...
            t_v, x_rg(1), x_rg(2), d_x, n_p, t_el, t_el / n_p);
end
fprintf('==================================================================================\n');