clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

lambda = 0.5;
L_list = 1:1:80;
num_L  = length(L_list);

funcs = {
    @(x) exp(-x.^2), ...
    @(x) 0.9 * exp(-x.^2 + x), ...
    @(x) sech(x), ...
    @(x) exp(-x.^2) .* sech(x)
};

func_names_console = {'Gaussian', 'Asymmetric Gaussian', 'sech', 'sech*Gaussian'};
func_legends_latex = {
    '$\mathrm{e}^{-x^2}$', ...
    '$0.9\,\mathrm{e}^{-x^2+x}$', ...
    '$\mathrm{sech}\,x$', ...
    '$\mathrm{e}^{-x^2}\mathrm{sech}\,x$'
};
num_funcs = length(funcs);

rho_all = zeros(num_funcs, num_L);

fprintf('=====================================================\n');
fprintf('  Truncation Convergence Test (lambda = %.2f)\n', lambda);
fprintf('=====================================================\n\n');

for f_idx = 1:num_funcs
    f = funcs{f_idx};
    fprintf('>>> Testing profile: %s (%d/%d)\n', func_names_console{f_idx}, f_idx, num_funcs);
    rho_vals = zeros(1, num_L);
    
    for i = 1:num_L
        L = L_list(i);
        
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
    rho_ref = rho_vals(end);
    err_rho = abs(rho_vals - rho_ref);
    
    fprintf('+------+---------------------+-------------+\n');
    fprintf('|  L   |      |rho|          |    Error    |\n');
    fprintf('+------+---------------------+-------------+\n');
    for idx = [5, 10, 20, 30, 40, 50, 60, 70, 80]
        L_val = L_list(idx);
        fprintf('| %4d | %17.10e | %10.2e |\n', L_val, abs(rho_vals(idx)), err_rho(idx));
    end
    fprintf('+------+---------------------+-------------+\n\n');
end

fprintf('=== Stability Assessment in Asymptotic Plateau (L >= 40) ===\n');
for f_idx = 1:num_funcs
    rho_stable = rho_all(f_idx, 40:end);
    range_rho = max(abs(rho_stable)) - min(abs(rho_stable));
    fprintf('%-24s: variation range = %.2e\n', func_names_console{f_idx}, range_rho);
end
fprintf('============================================================\n\n');

figure('Color', 'w', 'Position', [150, 150, 720, 520]);
hold on;

line_styles = {'-', '--', '-.', ':'};
colors      = {'#0072BD', '#D95319', '#2CA02C', '#111111'};
line_widths = [2.0, 1.8, 1.8, 2.0];

plot_idx = 1:(num_L - 1); 

for f_idx = 1:num_funcs
    err_this = abs(rho_all(f_idx, :) - rho_all(f_idx, end));
    plot(L_list(plot_idx), err_this(plot_idx), ...
         'LineStyle', line_styles{f_idx}, ...
         'Color', colors{f_idx}, ...
         'LineWidth', line_widths(f_idx));
end
hold off;

set(gca, 'YScale', 'log');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'LineWidth', 1.0);
box on; grid on;

xlabel('$L$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$E(L)$', 'Interpreter', 'latex', 'FontSize', 14);

legend(func_legends_latex, 'Interpreter', 'latex', 'FontSize', 12, ...
       'Location', 'southwest', 'Box', 'on');
xlim([1, 80]);
ylim([1e-13, 1e0]);