% =========================================================================
% Script: error_different_initial.m
% Description: Validation of the Chebyshev spectral collocation direct
%              scattering solver against exact closed-form expressions for
%              one-sided benchmark initial profiles.
%
% Generates:
%   - figure/different_initial_Serror.jpg    (Pointwise scatter plot)
%   - figure/different_initial_S_maxerror.jpg (Maximum absolute error bars)
%
% Reference:
%   "Numerical Inverse Scattering Transform for the Defocusing Nonlocal
%    Nonlinear Schrodinger Equation", Section 3.2, Figures 5-6.
% Repository: https://github.com/Bingtang-Orange/DNNLSE-NIST
% =========================================================================

clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% 1. Workspace Configuration and Directory Setup
% Resolve absolute path to guarantee reliable image export across platforms
script_path = mfilename('fullpath');
if ~isempty(script_path)
    output_dir = fullfile(fileparts(script_path), 'figure');
else
    output_dir = fullfile(pwd, 'figure');
end

if ~exist(output_dir, 'dir')
    [status, msg] = mkdir(output_dir);
    if ~status
        warning('Failed to create target directory: %s. Using current working directory.', msg);
        output_dir = pwd;
    end
end

%% 2. Spectral Grid and Computational Parameters
L = 40;                             % Domain truncation half-width
lambda_vec = -2.0:0.1:2.0;          % Continuous spectral grid (41 points)
num_lambda = length(lambda_vec);

%% 3. Benchmark Initial Profiles
% Profiles ordered by regularity at the origin: C0 -> C1 -> C2 -> Modulated C1
ICs = struct();

% Profile 1: C^0 regularity at x = 0
ICs(1).name       = '$s\mathrm{e}^{-s}\,(C^0)$';
ICs(1).name_plain = 's*exp(-s)';
ICs(1).f          = @(s) s .* exp(-s);
ICs(1).I1_expr    = @(lam) 1 ./ (1 - 2i*lam).^2;
ICs(1).I2_expr    = @(lam) 1 ./ (1 - 2i*lam).^2;

% Profile 2: C^1 regularity at x = 0
ICs(2).name       = '$s^2\mathrm{e}^{-s}\,(C^1)$';
ICs(2).name_plain = 's^2*exp(-s)';
ICs(2).f          = @(s) s.^2 .* exp(-s);
ICs(2).I1_expr    = @(lam) 2 ./ (1 - 2i*lam).^3;
ICs(2).I2_expr    = @(lam) 2 ./ (1 - 2i*lam).^3;

% Profile 3: C^2 regularity at x = 0
ICs(3).name       = '$s^3\mathrm{e}^{-s}\,(C^2)$';
ICs(3).name_plain = 's^3*exp(-s)';
ICs(3).f          = @(s) s.^3 .* exp(-s);
ICs(3).I1_expr    = @(lam) 6 ./ (1 - 2i*lam).^4;
ICs(3).I2_expr    = @(lam) 6 ./ (1 - 2i*lam).^4;

% Profile 4: Phase-modulated complex profile
ICs(4).name       = '$s^2\mathrm{e}^{-(1-\mathrm{i})s}$';
ICs(4).name_plain = 's^2*exp(-s)*exp(is)';
ICs(4).f          = @(s) s.^2 .* exp(-s) .* exp(1i*s);
ICs(4).I1_expr    = @(lam) 2 ./ (1 + 1i - 2i*lam).^3;
ICs(4).I2_expr    = @(lam) 2 ./ (1 - 1i - 2i*lam).^3;

num_ICs = length(ICs);

err_s11 = zeros(num_ICs, num_lambda);
err_s21 = zeros(num_ICs, num_lambda);

fprintf('=================================================================\n');
fprintf('  Direct Scattering Transform: Accuracy on One-Sided Benchmarks\n');
fprintf('=================================================================\n\n');

%% 4. Numerical Solution of the Spatial Lax System
for idx = 1:num_ICs
    f = ICs(idx).f;
    fprintf('>>> Computing scattering data for: %s (%d/%d)\n', ICs(idx).name_plain, idx, num_ICs);
    
    % Assemble piecewise smooth potentials on decoupled half-lines
    pot = struct();
    pot.L             = L;
    pot.q_left        = chebfun(0, [-L, 0]);
    pot.q_neg_c_left  = chebfun(@(x) conj(f(-x)), [-L, 0], 'splitting', 'on');
    pot.q_right       = chebfun(@(x) f(x), [0, L], 'splitting', 'on');
    pot.q_neg_c_right = chebfun(0, [0, L]);
    
    s11_arr = zeros(1, num_lambda);
    s21_arr = zeros(1, num_lambda);
    
    for k = 1:num_lambda
        lam = lambda_vec(k);
        S = compute_S(lam, pot);
        s11_arr(k) = S(1,1);
        s21_arr(k) = S(2,1);
    end
    
    % Evaluate closed-form reference scattering coefficients
    I1_ref = ICs(idx).I1_expr(lambda_vec);
    I2_ref = ICs(idx).I2_expr(lambda_vec);
    s11_ref = 1 + I1_ref .* I2_ref;
    s21_ref = I1_ref;
    
    % Compute absolute discrepancies
    err_s11(idx, :) = abs(s11_arr - s11_ref);
    err_s21(idx, :) = abs(s21_arr - s21_ref);
end

%% 5. Console Summary of Maximum Absolute Errors
fprintf('\n==================== Error Summary Table ========================\n');
fprintf('%-25s %-18s %-18s\n', 'Benchmark Profile', 'max |Err(s11)|', 'max |Err(s21)|');
fprintf('-----------------------------------------------------------------\n');
for idx = 1:num_ICs
    fprintf('%-25s %-18.4e %-18.4e\n', ICs(idx).name_plain, max(err_s11(idx,:)), max(err_s21(idx,:)));
end
fprintf('=================================================================\n\n');

%% 6. Figure 1: Pointwise Absolute Errors (Pure Scatter Plot)
fig_err = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
colors  = {'#0072BD', '#D95319', '#2CA02C', '#7E2F8E'};
markers = {'o', 's', '^', 'd'};

% Top subplot: Absolute error for s11(lambda)
subplot(2, 1, 1);
hold on;
for idx = 1:num_ICs
    scatter(lambda_vec, err_s11(idx, :), 36, ...
            'Marker', markers{idx}, ...
            'MarkerFaceColor', colors{idx}, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.85, ...
            'DisplayName', ICs(idx).name);
end
hold off;
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
box on; grid on;
ylabel('$\mathcal{E}(s_{11})$', 'Interpreter', 'latex', 'FontSize', 13);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');
xlim([-2, 2]); ylim([1e-12, 1e-8]);

% Bottom subplot: Absolute error for s21(lambda)
subplot(2, 1, 2);
hold on;
for idx = 1:num_ICs
    scatter(lambda_vec, err_s21(idx, :), 36, ...
            'Marker', markers{idx}, ...
            'MarkerFaceColor', colors{idx}, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.85, ...
            'DisplayName', ICs(idx).name);
end
hold off;
set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', 11);
box on; grid on;
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$\mathcal{E}(s_{21})$', 'Interpreter', 'latex', 'FontSize', 13);
legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast', 'Box', 'on');
xlim([-2, 2]); ylim([1e-12, 1e-8]);

% Export Figure 1 to the figure directory (300 DPI)
fig1_file = fullfile(output_dir, 'different_initial_Serror.jpg');
exportgraphics(fig_err, fig1_file, 'Resolution', 300);
fprintf('>> Exported: %s\n', fig1_file);

%% 7. Figure 2: Maximum Absolute Error Comparison (Grouped Bar Chart)
fig_bar = figure('Color', 'w', 'Position', [150, 150, 850, 420]);

max_errors = [max(err_s11, [], 2), max(err_s21, [], 2)];

b = bar(1:num_ICs, max_errors, 0.7);
b(1).FaceColor = [0.20, 0.40, 0.80];
b(2).FaceColor = [0.85, 0.35, 0.25];

set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:num_ICs, 'XTickLabel', {ICs.name}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('Maximum Absolute Error', 'Interpreter', 'latex', 'FontSize', 13);
legend({'$\max\,\mathcal{E}(s_{11})$', '$\max\,\mathcal{E}(s_{21})$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
box on; grid on;
ylim([1e-12, 1e-8]);

% Export Figure 2 to the figure directory (300 DPI)
fig2_file = fullfile(output_dir, 'different_initial_S_maxerror.jpg');
exportgraphics(fig_bar, fig2_file, 'Resolution', 300);
fprintf('>> Exported: %s\n', fig2_file);