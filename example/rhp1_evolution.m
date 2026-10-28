clear; clc; close all;

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
L = 5;
n = 400;

x_range_max = 30;
num_pts = 601;
x_grid = linspace(-x_range_max, x_range_max, num_pts).';
t_list = [0, 1, 2];

x_I = flipud(cos(pi*(0:n-1)'/(n-1)));
def = create_curve_def('interval', -L, L);
pts = curve_mapping(def, x_I);
lambda = pts{1};

q0 = chebfun(q0_handle, [-1.5*L, 1.5*L]);

fprintf('Computing scattering data (L = %.1f, n = %d)...\n', L, n);
tic;
[rho, rho_hat] = compute_reflection(lambda, q0, L);
fprintf('Done in %.2f s.\n\n', toc);

Q = cell(length(t_list), 1);
for kt = 1:length(t_list)
    t_try = t_list(kt);
    fprintf('Solving t = %.1f...\n', t_try);
    tic;
    qv = zeros(num_pts, 1);
    for j = 1:num_pts
        qv(j) = solve_rhp1(x_grid(j), t_try, lambda, rho, rho_hat, L, n);
    end
    Q{kt} = qv;
    fprintf('  Elapsed time: %.2f s\n', toc);
    fprintf('  max|q| = %.6e\n\n', max(abs(qv)));
end

figure('Position', [100, 100, 950, 350]);
plot_full(x_grid, Q{1}, 0, -20, 20, q0_handle);

figure('Position', [100, 100, 950, 1050]);

subplot(3, 1, 1);
plot_full(x_grid, Q{2}, 1, -20, 20);

subplot(3, 1, 2);
plot_local(x_grid, Q{2}, 1, -15, -5);

subplot(3, 1, 3);
plot_local(x_grid, Q{2}, 1, 5, 15);

sgtitle('$q(x, t=1)$', 'Interpreter', 'latex', 'FontSize', 13);

figure('Position', [100, 100, 950, 1050]);

subplot(3, 1, 1);
plot_full(x_grid, Q{3}, 2, -30, 30);

subplot(3, 1, 2);
plot_local(x_grid, Q{3}, 2, -25, -10);

subplot(3, 1, 3);
plot_local(x_grid, Q{3}, 2, 10, 25);

sgtitle('$q(x, t=2)$', 'Interpreter', 'latex', 'FontSize', 13);

function plot_full(x_grid, qv, t_val, x_lo, x_hi, q0_handle)
    mask = x_grid >= x_lo & x_grid <= x_hi;
    xx = x_grid(mask);
    qq = qv(mask);
    
    plot(xx, real(qq), 'r-',  'LineWidth', 1.4, 'DisplayName', 'Re q'); hold on;
    plot(xx, imag(qq), 'b--', 'LineWidth', 1.4, 'DisplayName', 'Im q');
    
    if nargin >= 6 && ~isempty(q0_handle) && t_val == 0
        dx_exact = 0.5;
        x_exact = (x_lo:dx_exact:x_hi).';
        q0_exact = q0_handle(x_exact);
        plot(x_exact, q0_exact, 'ko', 'MarkerSize', 4, 'LineWidth', 0.8, ...
             'MarkerFaceColor', 'none', 'DisplayName', 'Exact initial condition');
    end
    
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$t = %g$, full domain $x \\in [%d, %d]$', t_val, x_lo, x_hi), ...
          'Interpreter', 'latex', 'FontSize', 12);
    legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
    grid on; box on;
    xlim([x_lo, x_hi]);
end

function plot_local(x_grid, qv, t_val, x_lo, x_hi)
    mask = x_grid >= x_lo & x_grid <= x_hi;
    xx = x_grid(mask);
    qq = qv(mask);
    plot(xx, real(qq), 'r-',  'LineWidth', 1.4, 'DisplayName', 'Re q'); hold on;
    plot(xx, imag(qq), 'b--', 'LineWidth', 1.4, 'DisplayName', 'Im q');
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$q(x, t)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$t = %g$, local view $x \\in [%d, %d]$', t_val, x_lo, x_hi), ...
          'Interpreter', 'latex', 'FontSize', 12);
    legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'best');
    grid on; box on;
    xlim([x_lo, x_hi]);
    ymax = max(abs([real(qq); imag(qq)]));
    if ymax > 0
        ylim([-1.15*ymax, 1.15*ymax]);
    end
end