% RUN_FIG_RHP1_EVOLUTION
% Computes short-time spatio-temporal solutions of DNNLSE via RHP I solver.
% Benchmark initial profile: q0(x) = 0.9 * exp(-x^2 + x) at t = 0, 1, 2.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

fprintf('===================================================================\n');
fprintf(' RHP I Direct Solver: Short-Time Dynamics for q0(x) = 0.9 exp(-x^2+x)\n');
fprintf('===================================================================\n\n');

% ---------- Simulation Parameters ----------
q0_handle = @(x) 0.9 * exp(-x.^2 + x);

L_spec    = 5.0;     % Spectral truncation half-length for RHP I jump contour
L_spatial = 40.0;    % Physical domain truncation length for direct scattering
n_nodes   = 400;     % Chebyshev collocation points on spectral contour

% Spatial evaluation domain
x_max   = 30.0;
num_pts = 601;
x_grid  = linspace(-x_max, x_max, num_pts).';

% Temporal snapshots
t_list = [0.0, 1.0, 2.0];

% ---------- Setup Spectral Grid and Precompute Scattering Data ----------
x_I = flipud(cos(pi * (0:n_nodes-1)' / (n_nodes - 1)));
def = create_curve_def('interval', -L_spec, L_spec);
pts = curve_mapping(def, x_I);
lambda = pts{1};

fprintf('Computing direct scattering data (L_spatial = %.1f, n = %d) ...\n', ...
        L_spatial, n_nodes);
tic;
[rho, rho_hat] = compute_reflection(lambda, q0_handle, L_spatial);
fprintf('Scattering data generated in %.2f s.\n\n', toc);

% ---------- Solve RHP I Across Temporal Snapshots ----------
Q_solutions = cell(length(t_list), 1);

for kt = 1:length(t_list)
    t_val = t_list(kt);
    fprintf('--> Solving RHP I at snapshot t = %.1f ...\n', t_val);
    tic;
    qv = zeros(num_pts, 1);
    for j = 1:num_pts
        qv(j) = solve_rhp1(x_grid(j), t_val, lambda, rho, rho_hat, L_spec, n_nodes);
    end
    Q_solutions{kt} = qv;
    fprintf('    Solved %d points in %.2f s (max|q| = %.4e).\n', ...
            num_pts, toc, max(abs(qv)));
end

% =========================================================================
% Figure 1: Initial Profile at t = 0 (Single Panel)
% =========================================================================
figure('Color', 'w', 'Position', [100, 100, 850, 320]);
plot_full(x_grid, Q_solutions{1}, 0, -20, 20);

% =========================================================================
% Figure 2: Snapshot at t = 1 (Full Domain and Local Wavepacket Zoom)
% =========================================================================
figure('Color', 'w', 'Position', [120, 80, 850, 750]);

subplot(3, 1, 1);
plot_full(x_grid, Q_solutions{2}, 1, -20, 20);

subplot(3, 1, 2);
plot_local(x_grid, Q_solutions{2}, 1, -15, -5);

subplot(3, 1, 3);
plot_local(x_grid, Q_solutions{2}, 1, 5, 15);

% =========================================================================
% Figure 3: Snapshot at t = 2 (Full Domain and Dispersive Tails Zoom)
% =========================================================================
figure('Color', 'w', 'Position', [140, 60, 850, 750]);

subplot(3, 1, 1);
plot_full(x_grid, Q_solutions{3}, 2, -30, 30);

subplot(3, 1, 2);
plot_local(x_grid, Q_solutions{3}, 2, -25, -10);

subplot(3, 1, 3);
plot_local(x_grid, Q_solutions{3}, 2, 10, 25);

fprintf('\nAll dynamic profiles successfully computed and visualized.\n');

% =========================================================================
% Helper Plotting Functions
% =========================================================================
function plot_full(x_grid, qv, t_val, x_lo, x_hi)
    mask = (x_grid >= x_lo) & (x_grid <= x_hi);
    xx = x_grid(mask);
    qq = qv(mask);

    plot(xx, real(qq), 'r-',  'LineWidth', 1.4); hold on;
    plot(xx, imag(qq), 'b--', 'LineWidth', 1.4); hold off;
    grid on; box on;
    set(gca, 'FontSize', 10, 'LineWidth', 1.0);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$t = %.0f$, full domain $x \\in [%d, %d]$', t_val, x_lo, x_hi), ...
          'Interpreter', 'latex', 'FontSize', 11);
    legend({'$\mathrm{Re}\,q$', '$\mathrm{Im}\,q$'}, ...
           'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    xlim([x_lo, x_hi]);
end

function plot_local(x_grid, qv, t_val, x_lo, x_hi)
    mask = (x_grid >= x_lo) & (x_grid <= x_hi);
    xx = x_grid(mask);
    qq = qv(mask);

    plot(xx, real(qq), 'r-',  'LineWidth', 1.4); hold on;
    plot(xx, imag(qq), 'b--', 'LineWidth', 1.4); hold off;
    grid on; box on;
    set(gca, 'FontSize', 10, 'LineWidth', 1.0);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$t = %.0f$, local zoom $x \\in [%d, %d]$', t_val, x_lo, x_hi), ...
          'Interpreter', 'latex', 'FontSize', 11);
    legend({'$\mathrm{Re}\,q$', '$\mathrm{Im}\,q$'}, ...
           'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);
    xlim([x_lo, x_hi]);

    ymax = max(abs([real(qq); imag(qq)]));
    if ymax > 0
        ylim([-1.15 * ymax, 1.15 * ymax]);
    end
end