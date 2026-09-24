% =========================================================================
% Script: plot_reflection_coefficients.m
% Description: Computes and plots the spectral reflection coefficients
%              rho(lambda) and rho_hat(lambda) for the benchmark Gaussian
%              initial condition q(x) = exp(-x^2) across the real continuous
%              spectrum lambda in [-2, 2].
%
% Generates:
%   - figure/reflection_coeffcient_lambda.jpg (Figure 7 in manuscript)
%
% Reference:
%   "Numerical Inverse Scattering Transform for the Defocusing Nonlocal
%    Nonlinear Schrodinger Equation", Section 3.2, Eq. (8) and Figure 7.
% Repository: https://github.com/Bingtang-Orange/DNNLSE-NIST
% =========================================================================

clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% 1. Workspace Configuration and Output Directory
script_path = mfilename('fullpath');
if ~isempty(script_path)
    output_dir = fullfile(fileparts(script_path), 'figure');
else
    output_dir = fullfile(pwd, 'figure');
end

if ~exist(output_dir, 'dir')
    [status, msg] = mkdir(output_dir);
    if ~status
        warning('Failed to create target directory: %s. Using current folder.', msg);
        output_dir = pwd;
    end
end

%% 2. Computational Setup and Parameters
L = 20;                             % Domain truncation half-width
lambda_vec = -2.0:0.05:2.0;         % Real spectral grid (81 points)
num_lambda = length(lambda_vec);

% Gaussian benchmark profile q(x) = exp(-x^2)
q_func = @(x) exp(-x.^2);

% Assemble piecewise potentials on decoupled half-lines
pot = struct();
pot.L             = L;
pot.q_left        = chebfun(@(x) q_func(x),        [-L, 0]);
pot.q_neg_c_left  = chebfun(@(x) conj(q_func(-x)), [-L, 0]);
pot.q_right       = chebfun(@(x) q_func(x),        [0, L]);
pot.q_neg_c_right = chebfun(@(x) conj(q_func(-x)), [0, L]);

%% 3. Numerical Evaluation of the Direct Scattering Transform
s11_arr     = zeros(1, num_lambda);
s12_arr     = zeros(1, num_lambda);
s21_arr     = zeros(1, num_lambda);
s22_arr     = zeros(1, num_lambda);
rho_arr     = zeros(1, num_lambda);
rho_hat_arr = zeros(1, num_lambda);

fprintf('=================================================================\n');
fprintf('  Direct Scattering Evaluation: Gaussian Benchmark Profile\n');
fprintf('=================================================================\n');
fprintf('Computing scattering data across lambda in [%.2f, %.2f]...\n', ...
        lambda_vec(1), lambda_vec(end));

for k = 1:num_lambda
    lam = lambda_vec(k);
    S = compute_S(lam, pot);
    
    s11_arr(k) = S(1,1);
    s12_arr(k) = S(1,2);
    s21_arr(k) = S(2,1);
    s22_arr(k) = S(2,2);
    
    % Reflection coefficients strictly according to Eq. (8) in manuscript
    rho_arr(k)     = s12_arr(k) / s22_arr(k);
    rho_hat_arr(k) = s21_arr(k) / s11_arr(k);
end

%% 4. Figure: Reflection Coefficients (Fixed LaTeX Syntax)
fig = figure('Color', 'w', 'Position', [150, 100, 800, 600]);

% Color and line specifications
c_real = '#0072BD';   % Deep blue
c_imag = '#D95319';   % Brick red
c_mod  = '#111111';   % Solid dark

% Subplot 1: Left-incident reflection coefficient rho(lambda)
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
xlim([-2, 2]);
ylabel('$\rho(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');

% Subplot 2: Right-incident reflection coefficient rho_hat(lambda)
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
xlim([-2, 2]);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{\rho}(\lambda)$', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');

% Export Figure 7 to the figure directory (300 DPI)
export_path = fullfile(output_dir, 'reflection_coeffcient_lambda.jpg');
exportgraphics(fig, export_path, 'Resolution', 300);
fprintf('>> Exported successfully to: %s\n', export_path);

%% 5. Console Summary and Theoretical PT-Symmetry Check
sym_error = max(abs(rho_arr - conj(rho_hat_arr)));

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

fprintf('PT-Symmetry Verification (rho == conj(rho_hat)):\n');
fprintf('  Max discrepancy max |rho - conj(rho_hat)|: %.4e\n', sym_error);
fprintf('=================================================================\n');