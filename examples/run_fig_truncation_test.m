% RUN_FIG_TRUNCATION_TEST 
% Evaluates reflection coefficient convergence against truncation length L.
% Reproduces Figure for Section 3.2.2.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

% ---- Simulation Parameters ----
lambda = 0.5;
L_list = 1:1:80;
num_L  = length(L_list);

% Benchmark potentials
funcs = {
    @(x) exp(-x.^2),                    % Gaussian
    @(x) exp(-x.^2 + 1i*x),             % Phase-shifted Gaussian
    @(x) sech(x),                       % Hyperbolic secant
    @(x) exp(-x.^2) .* sech(x)          % Modulated sech
};
func_names = {'Gaussian', 'Phase-shifted Gaussian', 'sech', 'sech*Gaussian'};
num_funcs  = length(funcs);

rho_all = zeros(num_funcs, num_L);

fprintf('=====================================================\n');
fprintf(' Truncation Convergence Test at lambda = %.2f\n', lambda);
fprintf('=====================================================\n\n');

% ---- Computation Loop ----
for f_idx = 1:num_funcs
    f = funcs{f_idx};
    fprintf('--> Computing for: %s ...\n', func_names{f_idx});
    rho_vals = zeros(1, num_L);
    
    for i = 1:num_L
        L = L_list(i);
        % Use the high-level direct scattering interface
        [rho_vals(i), ~] = compute_reflection(lambda, f, L);
    end
    
    rho_all(f_idx, :) = rho_vals;
    rho_ref = rho_vals(end); % Reference at L = 80
    err_rho = abs(rho_vals - rho_ref);
    
    % Display verification table
    fprintf('+------+---------------------+-------------+\n');
    fprintf('|  L   |      |rho|          |   Error     |\n');
    fprintf('+------+---------------------+-------------+\n');
    sample_indices = [5, 10, 20, 30, 40, 50, 60, 70, 80];
    for idx = sample_indices
        L = L_list(idx);
        fprintf('| %4d | %17.10e | %10.2e |\n', L, abs(rho_vals(idx)), err_rho(idx));
    end
    fprintf('+------+---------------------+-------------+\n\n');
end

% ---- Plotting: Error vs. Domain Length L (Paper Figure) ----
figure('Color', 'w', 'Position', [150, 150, 700, 480]);
hold on;
line_styles = {'-', '--', '-.', ':'};
colors      = {'b', 'r', [0 0.5 0], 'k'}; 

for f_idx = 1:num_funcs
    err_this = abs(rho_all(f_idx, :) - rho_all(f_idx, end));
    plot(L_list, err_this, 'LineWidth', 1.8, ...
         'LineStyle', line_styles{f_idx}, 'Color', colors{f_idx});
end
hold off;

set(gca, 'YScale', 'log', 'FontSize', 11, 'LineWidth', 1.1);
xlabel('Truncation length $L$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$|\rho(L) - \rho(80)|$', 'Interpreter', 'latex', 'FontSize', 13);
title('Convergence with Respect to Domain Truncation $L$', 'Interpreter', 'latex', 'FontSize', 13);
legend(func_names, 'Location', 'northeast', 'Interpreter', 'none', 'FontSize', 10);
xlim([0, 80]);
ylim([1e-13, 1e0]);
grid on; box on;