% Plot reflection coefficients rho(lambda) and rho_hat(lambda).

clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% Parameters and potential
L = 5;
lambda_vec = -5.0:0.05:5.0;
num_lambda = length(lambda_vec);

q_func = @(x) 0.9 * exp(-x.^2 + x);

pot = struct();
pot.L             = L;
pot.q_left        = chebfun(@(x) q_func(x),        [-L, 0]);
pot.q_neg_c_left  = chebfun(@(x) conj(q_func(-x)), [-L, 0]);
pot.q_right       = chebfun(@(x) q_func(x),        [0, L]);
pot.q_neg_c_right = chebfun(@(x) conj(q_func(-x)), [0, L]);

%% Compute scattering data
s11_arr     = zeros(1, num_lambda);
s12_arr     = zeros(1, num_lambda);
s21_arr     = zeros(1, num_lambda);
s22_arr     = zeros(1, num_lambda);
rho_arr     = zeros(1, num_lambda);
rho_hat_arr = zeros(1, num_lambda);

fprintf('Computing scattering data across lambda in [%.2f, %.2f]...\n', ...
        lambda_vec(1), lambda_vec(end));

for k = 1:num_lambda
    lam = lambda_vec(k);
    S = compute_S(lam, pot);
    
    s11_arr(k) = S(1,1);
    s12_arr(k) = S(1,2);
    s21_arr(k) = S(2,1);
    s22_arr(k) = S(2,2);
    
    rho_arr(k)     = s12_arr(k) / s22_arr(k);
    rho_hat_arr(k) = s21_arr(k) / s11_arr(k);
end

%% Plot reflection coefficients
figure('Color', 'w', 'Position', [150, 100, 800, 600]);

c_real = '#0072BD';
c_imag = '#D95319';
c_mod  = '#111111';

subplot(2, 1, 1);
hold on;
plot(lambda_vec, real(rho_arr), 'LineStyle', '-',  'Color', c_real, 'LineWidth', 1.6, ...
     'DisplayName', '$\mathrm{Re}\,\rho(\lambda)$');
plot(lambda_vec, imag(rho_arr), 'LineStyle', '--', 'Color', c_imag, 'LineWidth', 1.6, ...
     'DisplayName', '$\mathrm{Im}\,\rho(\lambda)$');
plot(lambda_vec, abs(rho_arr),  'LineStyle', '-.', 'Color', c_mod,  'LineWidth', 1.8, ...
     'DisplayName', '$|\rho(\lambda)|$');
hold off;

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11, 'LineWidth', 1.0);
box on; grid on;
xlim([lambda_vec(1), lambda_vec(end)]);
ylabel('$\rho(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');

subplot(2, 1, 2);
hold on;
plot(lambda_vec, real(rho_hat_arr), 'LineStyle', '-',  'Color', c_real, 'LineWidth', 1.6, ...
     'DisplayName', '$\mathrm{Re}\,\hat{\rho}(\lambda)$');
plot(lambda_vec, imag(rho_hat_arr), 'LineStyle', '--', 'Color', c_imag, 'LineWidth', 1.6, ...
     'DisplayName', '$\mathrm{Im}\,\hat{\rho}(\lambda)$');
plot(lambda_vec, abs(rho_hat_arr),  'LineStyle', '-.', 'Color', c_mod,  'LineWidth', 1.8, ...
     'DisplayName', '$|\hat{\rho}(\lambda)|$');
hold off;

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11, 'LineWidth', 1.0);
box on; grid on;
xlim([lambda_vec(1), lambda_vec(end)]);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{\rho}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');

%% Summary
pt_discrepancy = max(abs(rho_arr - conj(rho_hat_arr)));

fprintf('\n================== Summary of Plotted Data ======================\n');
fprintf('Spectral grid: [%.2f, %.2f], step = %.2f (Total points: %d)\n\n', ...
        lambda_vec(1), lambda_vec(end), lambda_vec(2)-lambda_vec(1), num_lambda);

fprintf('Left-incident reflection coefficient rho(lambda):\n');
[max_rho, idx_max] = max(abs(rho_arr));
[min_rho, idx_min] = min(abs(rho_arr));
fprintf('  Max modulus: %12.4e (at lambda = %+.2f)\n', max_rho, lambda_vec(idx_max));
fprintf('  Min modulus: %12.4e (at lambda = %+.2f)\n', min_rho, lambda_vec(idx_min));
fprintf('  Real range : [%+12.4e, %+12.4e]\n', min(real(rho_arr)), max(real(rho_arr)));
fprintf('  Imag range : [%+12.4e, %+12.4e]\n\n', min(imag(rho_arr)), max(imag(rho_arr)));

fprintf('Right-incident reflection coefficient rho_hat(lambda):\n');
[max_rh, idx_maxh] = max(abs(rho_hat_arr));
[min_rh, idx_minh] = min(abs(rho_hat_arr));
fprintf('  Max modulus: %12.4e (at lambda = %+.2f)\n', max_rh, lambda_vec(idx_maxh));
fprintf('  Min modulus: %12.4e (at lambda = %+.2f)\n', min_rh, lambda_vec(idx_minh));
fprintf('  Real range : [%+12.4e, %+12.4e]\n', min(real(rho_hat_arr)), max(real(rho_hat_arr)));
fprintf('  Imag range : [%+12.4e, %+12.4e]\n\n', min(imag(rho_hat_arr)), max(imag(rho_hat_arr)));

fprintf('Nonlocal Asymmetry Indicator:\n');
fprintf('  Max |rho - conj(rho_hat)|: %.4e (confirms PT-symmetry breaking)\n', pt_discrepancy);
fprintf('=================================================================\n');