clear; clc; close all;

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

ratio_c = -3.0; 
lambda0_const = ratio_c / 4.0;

[r0, rh0] = compute_reflection(lambda0_const, q0, 8.0);
det_fac0 = 1.0 - r0 * rh0;

nu0 = -(1 / (2 * pi)) * log(det_fac0);
nu_R = real(nu0);
nu_I = imag(nu0);

alpha_theory = -0.50 - nu_I;

fprintf('====================================================================\n');
fprintf('  Asymptotic Decay Verification: x = %g * t (lambda0 = %g)\n', ratio_c, lambda0_const);
fprintf('====================================================================\n');
fprintf('  1 - rho*rho_hat  = %.6f + %.6fi\n', real(det_fac0), imag(det_fac0));
fprintf('  arg(1 - rho*rho) = %.6f rad\n', angle(det_fac0));
fprintf('  nu_R             = %.6f\n', nu_R);
fprintf('  nu_I             = %.6f\n', nu_I);
fprintf('  Theoretical decay slope alpha_theory = %.6f\n\n', alpha_theory);

num_t = 41;
t_grid = logspace(2, 4, num_t).';
x_grid = ratio_c * t_grid;

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

q_vals = complex(zeros(num_t, 1));
zs_records = zeros(num_t, 1);

fprintf('Solving RHP IV for t in [100, 10000]...\n');
t_start = tic;
parfor j = 1:num_t
    tj = t_grid(j);
    xj = x_grid(j);
    
    [qv, info] = solve_rhp4(xj, tj, q0, 8.0, 8.0, 24, 20);
    q_vals(j) = qv;
    zs_records(j) = info.zero_sum_max;
end
t_cost = toc(t_start);

fprintf('Computation completed in %.2f s (%.2f s/point).\n', t_cost, t_cost / num_t);
fprintf('Maximum zero-sum residual: %.2e\n\n', max(zs_records));

abs_q = abs(q_vals);
log_t = log10(t_grid);
log_q = log10(abs_q);

idx_asymp = t_grid >= 300;
fit_poly = polyfit(log_t(idx_asymp), log_q(idx_asymp), 1);
fit_slope = fit_poly(1);
fitted_log_q = polyval(fit_poly, log_t);

ref_log_q = alpha_theory * (log_t - log_t(end)) + log_q(end);

fprintf('================== Fitting vs Theory ====================\n');
fprintf('  Numerical slope   : alpha_num    = %.6f\n', fit_slope);
fprintf('  Theoretical slope : alpha_theory = %.6f\n', alpha_theory);
fprintf('  Absolute error    : |diff|       = %.4e\n', abs(fit_slope - alpha_theory));
fprintf('=========================================================\n\n');

figure('Color', 'w', 'Position', [100, 150, 1050, 420]);

subplot(1, 2, 1);
plot(t_grid, abs_q, 'b.-', 'LineWidth', 1.3, 'MarkerSize', 9, ...
     'DisplayName', '$|q(-3t, t)|$');
grid on; box on;
xlim([min(t_grid), max(t_grid)]);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$|q(x, t)|$', 'Interpreter', 'latex', 'FontSize', 12);
title('(a) Solution Amplitude Decay along $x = -3t$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);

subplot(1, 2, 2);
plot(log_t, log_q, 'ro', 'MarkerSize', 5, 'LineWidth', 1.2, ...
     'DisplayName', 'Numerical data'); hold on;
plot(log_t, fitted_log_q, 'k-', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('Fitted Line (Slope = %.4f)', fit_slope));
plot(log_t, ref_log_q, 'b--', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('Theoretical (Slope = %.4f)', alpha_theory));

grid on; box on;
xlim([log_t(1), log_t(end)]);
ylim([min(log_q) - 0.08, max(log_q) + 0.08]);
xlabel('$\log_{10}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\log_{10}|q(x, t)|$', 'Interpreter', 'latex', 'FontSize', 12);
title('(b) Verification of Modulated Power Law', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);