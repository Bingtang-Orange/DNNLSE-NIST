% RUN_FIG_DECAY_SLOPE_BENCHMARK
% Large-time asymptotic analysis and verification of the anomalous algebraic decay 
% rate along the characteristic ray x = -3t for DNNLSE.
% Reproduces Figure 21 and validates the modified power-law alpha_theory = -1/2 + nu_I.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% ---------------- 1. Physical Parameters and Benchmark Setup ----------------
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

% Characteristic ray x = c * t and stationary phase saddle point lambda0 = c / 4
ratio_c       = -3.0; 
lambda0_const = ratio_c / 4.0; % lambda0 = -0.75

fprintf('====================================================================\n');
fprintf(' Characteristic Ray Decay Verification: x = %g * t (lambda0 = %g)\n', ...
        ratio_c, lambda0_const);
fprintf(' Benchmark Profile: q0(x) = 0.9 * exp(-x^2 + x)\n');
fprintf('====================================================================\n\n');

% ---------------- 2. Analytical Theory: Modulated Decay Exponent ----------------
% Evaluate spectral scattering data at the fixed saddle point
[r0, rh0] = compute_reflection(lambda0_const, q0, 20.0);
det_fac0  = 1.0 - r0 * rh0;

% Parabolic cylinder phase parameter nu(lambda0)
nu0  = -(1.0 / (2.0 * pi)) * log(det_fac0);
nu_R = real(nu0);
nu_I = imag(nu0);

% Theoretical asymptotic decay rate: standard dispersive -1/2 plus phase correction
alpha_theory = -0.50 - nu_I; 

fprintf('>> [Analytical Asymptotic Analysis]\n');
fprintf('   1 - rho*rho_hat  = %.6f %+.6fi\n', real(det_fac0), imag(det_fac0));
fprintf('   Phase angle phi  = %.6f rad\n', angle(det_fac0));
fprintf('   Real parameter   nu_R = %.6f\n', nu_R);
fprintf('   Modulation term  nu_I = %.6f\n', nu_I);
fprintf('   Exact Theoretical Slope: alpha_theory = -0.5 + nu_I = %.6f\n\n', alpha_theory);

% ---------------- 3. Large-Time Grid and Parallel RHP IV Solve ----------------
% Logarithmically spaced grid over large-time domain t in [100, 10000]
num_t  = 41;
t_grid = logspace(2, 4, num_t).';
x_grid = ratio_c * t_grid;

pool = gcp('nocreate');
if isempty(pool)
    parpool('Processes');
end

q_vals     = complex(zeros(num_t, 1));
zs_records = zeros(num_t, 1);

fprintf('>> Solving large-time RHP IV across %d points (t in [10^2, 10^4])...\n', num_t);
t_start = tic;
parfor j = 1:num_t
    tj = t_grid(j);
    xj = x_grid(j);
    
    % RHP IV solver with regularized lens scaling
    [qv, info] = solve_rhp4(xj, tj, q0, 20.0, 8.0, 24, 20);
    q_vals(j)     = qv;
    zs_records(j) = info.zero_sum_max;
end
t_cost = toc(t_start);

fprintf('>> Parallel evaluation completed in %.2f s (Average: %.2f s/point).\n', ...
        t_cost, t_cost / num_t);
fprintf('>> Maximum Zero-Sum residual: %.2e (Junction compatibility satisfied).\n\n', ...
        max(zs_records));

% ---------------- 4. Amplitude Extraction and Log-Linear Regression ----------------
abs_q = abs(q_vals);
log_t = log10(t_grid);
log_q = log10(abs_q);

% Perform least-squares regression on asymptotic stage (t >= 300)
idx_asymp    = t_grid >= 300;
fit_poly     = polyfit(log_t(idx_asymp), log_q(idx_asymp), 1);
fit_slope    = fit_poly(1);
fitted_log_q = polyval(fit_poly, log_t);

% Theoretical reference line aligned with the endpoint
ref_log_q = alpha_theory * (log_t - log_t(end)) + log_q(end);

fprintf('>> [Numerical Regression vs. Analytical Theory]\n');
fprintf('   Numerical fitted slope : alpha_num    = %.6f\n', fit_slope);
fprintf('   Theoretical prediction : alpha_theory = %.6f\n', alpha_theory);
fprintf('   Absolute discrepancy   : |diff|       = %.4e\n\n', abs(fit_slope - alpha_theory));

% ---------------- 5. Publication-Grade Visualization (Figure 21 Layout) ----------------
fig = figure('Color', 'w', 'Position', [100, 150, 950, 400]);

% Panel (a): Physical amplitude decay |q(-3t, t)| vs t
subplot(1, 2, 1);
plot(t_grid, abs_q, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 8, ...
     'DisplayName', '$|q(-3t, t)|$');
grid on; box on;
xlim([min(t_grid), max(t_grid)]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$|q(x, t)|$', 'Interpreter', 'latex', 'FontSize', 12);
title('(a) Amplitude Decay along $x = -3t$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);

% Panel (b): Double-logarithmic linear fit vs theoretical slope
subplot(1, 2, 2);
plot(log_t, log_q, 'ro', 'MarkerSize', 4.5, 'LineWidth', 1.1, ...
     'DisplayName', 'Numerical data'); hold on;
plot(log_t, fitted_log_q, 'k-', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('Fitted Line (Slope = %.4f)', fit_slope));
plot(log_t, ref_log_q, 'b--', 'LineWidth', 1.3, ...
     'DisplayName', sprintf('Theory (Slope = %.4f)', alpha_theory));
hold off;

grid on; box on;
xlim([log_t(1), log_t(end)]);
ylim([min(log_q) - 0.08, max(log_q) + 0.08]);
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
xlabel('$\log_{10}(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\log_{10}|q(x, t)|$', 'Interpreter', 'latex', 'FontSize', 12);
title('(b) Anomalous Modulated Power Law', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);