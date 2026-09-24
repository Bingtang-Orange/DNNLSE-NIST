% =========================================================================
% Script: convergence_test.m
% Description: Convergence analysis of the reflection coefficient with
%              respect to the domain truncation half-width L (from 1 to 80)
%              for four representative initial profiles.
%
% Generates:
%   - figure/error_of_reflection.jpg (Figure 4 in manuscript)
%
% Reference:
%   "Numerical Inverse Scattering Transform for the Defocusing Nonlocal
%    Nonlinear Schrodinger Equation", Section 3.2, Eq. (22) and Figure 4.
% Repository: https://github.com/Bingtang-Orange/DNNLSE-NIST
% =========================================================================

clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% 1. Workspace Configuration and Output Directory Setup
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

%% 2. Computational Parameters and Benchmark Functions
lambda = 0.5;                       % Fixed evaluation spectral parameter
L_list = 1:1:80;                    % Truncation lengths from 1 to 80
num_L  = length(L_list);

% Define four benchmark profiles matching the manuscript
funcs = {
    @(x) exp(-x.^2),                    % Gaussian
    @(x) exp(-x.^2 + 1i*x),             % Phase-shifted Gaussian
    @(x) sech(x),                       % Hyperbolic secant
    @(x) exp(-x.^2) .* sech(x)          % Gaussian-modulated sech
};

func_names_console = {'Gaussian', 'Phase-shifted Gaussian', 'sech', 'sech*Gaussian'};
func_legends_latex = {
    '$\mathrm{e}^{-x^2}$', ...
    '$\mathrm{e}^{-x^2+\mathrm{i}x}$', ...
    '$\mathrm{sech}\,x$', ...
    '$\mathrm{e}^{-x^2}\mathrm{sech}\,x$'
};
num_funcs = length(funcs);

% Preallocate reflection coefficient matrix
rho_all = zeros(num_funcs, num_L);

fprintf('=====================================================\n');
fprintf('  Truncation Convergence Test (lambda = %.2f)\n', lambda);
fprintf('=====================================================\n\n');

%% 3. Main Computation Loop across Profiles and Truncation Lengths
for f_idx = 1:num_funcs
    f = funcs{f_idx};
    fprintf('>>> Testing profile: %s (%d/%d)\n', func_names_console{f_idx}, f_idx, num_funcs);
    rho_vals = zeros(1, num_L);
    
    for i = 1:num_L
        L = L_list(i);
        
        % Assemble piecewise Chebfun potentials on decoupled half-lines
        pot = struct();
        pot.L             = L;
        pot.q_left        = chebfun(@(x) f(x),        [-L, 0], 'splitting', 'on', 'maxLength', 3000);
        pot.q_neg_c_left  = chebfun(@(x) conj(f(-x)), [-L, 0], 'splitting', 'on', 'maxLength', 3000);
        pot.q_right       = chebfun(@(x) f(x),        [0, L],  'splitting', 'on', 'maxLength', 3000);
        pot.q_neg_c_right = chebfun(@(x) conj(f(-x)), [0, L],  'splitting', 'on', 'maxLength', 3000);
        
        S = compute_S(lambda, pot);
        rho_vals(i) = S(1,2) / S(2,2);
    end
    
    rho_all(f_idx, :) = rho_vals;
    rho_ref = rho_vals(end); % Numerical reference at L = 80
    err_rho = abs(rho_vals - rho_ref);
    
    % Print summary table for representative domain sizes
    fprintf('+------+---------------------+-------------+\n');
    fprintf('|  L   |      |rho|          |    Error    |\n');
    fprintf('+------+---------------------+-------------+\n');
    for idx = [5, 10, 20, 30, 40, 50, 60, 70, 80]
        L_val = L_list(idx);
        fprintf('| %4d | %17.10e | %10.2e |\n', L_val, abs(rho_vals(idx)), err_rho(idx));
    end
    fprintf('+------+---------------------+-------------+\n\n');
end

%% 4. Stability Assessment in the Plateau Regime (L >= 40)
fprintf('=== Stability Assessment in Asymptotic Plateau (L >= 40) ===\n');
for f_idx = 1:num_funcs
    rho_stable = rho_all(f_idx, 40:end);
    range_rho = max(abs(rho_stable)) - min(abs(rho_stable));
    fprintf('%-24s: variation range = %.2e\n', func_names_console{f_idx}, range_rho);
end
fprintf('============================================================\n\n');

%% 5. Figure 1: Truncation Error (Strictly Preserving Line Styles & Colors)
fig1 = figure('Color', 'w', 'Position', [150, 150, 720, 520]);
hold on;

% Line appearance settings preserved exactly as configured
line_styles = {'-', '--', '-.', ':'};
colors      = {'#0072BD', '#D95319', '#2CA02C', '#111111'};
line_widths = [2.0, 1.8, 1.8, 2.0];

% Plot up to L = 79 to prevent log-scale singularity at L = 80 (Error = 0)
plot_idx = 1:(num_L - 1); 

for f_idx = 1:num_funcs
    err_this = abs(rho_all(f_idx, :) - rho_all(f_idx, end));
    plot(L_list(plot_idx), err_this(plot_idx), ...
         'LineStyle', line_styles{f_idx}, ...
         'Color', colors{f_idx}, ...
         'LineWidth', line_widths(f_idx));
end
hold off;

% Publication-quality axes configuration
set(gca, 'YScale', 'log');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'LineWidth', 1.0);
box on; grid on;

xlabel('$L$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E(L) = |\rho(L) - \rho(80)|$', 'Interpreter', 'latex', 'FontSize', 14);

legend(func_legends_latex, 'Interpreter', 'latex', 'FontSize', 12, ...
       'Location', 'southwest', 'Box', 'on');
xlim([1, 80]);
ylim([1e-13, 1e0]);

% Export Figure 4 to the figure directory (300 DPI)
fig1_file = fullfile(output_dir, 'error_of_reflection.jpg');
exportgraphics(fig1, fig1_file, 'Resolution', 300);
fprintf('>> Figure successfully exported to: %s\n', fig1_file);