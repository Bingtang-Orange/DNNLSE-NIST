% =========================================================================
% Script: generate_reflection_table.m
% Description: Validates the direct scattering solver on the benchmark profile
%              q_0(x) = x^2 * exp(-(1-i)*x) * Theta(x) against exact closed-form
%              expressions, generating the LaTeX source for Table 1.
%
% Output:
%   - Console: Formatted error statistics (Max, Min, Mean)
%   - File:    reflection_table.tex (Publication-ready LaTeX table)
%
% Reference:
%   "Numerical Inverse Scattering Transform for the Defocusing Nonlocal
%    Nonlinear Schrodinger Equation", Section 3.2, Table 1.
% Repository: https://github.com/Bingtang-Orange/DNNLSE-NIST
% =========================================================================

clear; clc;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% 1. Computational Setup and Parameters
L = 40;                             % Domain truncation half-width
lambda_vec = -2.0:0.5:2.0;          % Evaluation spectral nodes
num_lambda = length(lambda_vec);

% Benchmark profile: f(s) = s^2 * exp(-s) * exp(i*s) = s^2 * exp(-(1-i)*s)
f = @(s) s.^2 .* exp(-s) .* exp(1i*s);

% Assemble piecewise potentials on decoupled half-lines
pot = struct();
pot.L             = L;
pot.q_left        = chebfun(0, [-L, 0]);
pot.q_neg_c_left  = chebfun(@(x) conj(f(-x)), [-L, 0], 'splitting', 'on');
pot.q_right       = chebfun(@(x) f(x), [0, L], 'splitting', 'on');
pot.q_neg_c_right = chebfun(0, [0, L]);

%% 2. Direct Scattering Computation via Chebfun
s11_num = zeros(1, num_lambda);
s21_num = zeros(1, num_lambda);

fprintf('Computing numerical scattering data via Chebyshev collocation...\n');
for k = 1:num_lambda
    lam = lambda_vec(k);
    S = compute_S(lam, pot);
    s11_num(k) = S(1,1);
    s21_num(k) = S(2,1);
end

%% 3. Closed-Form Exact Scattering Data
% Closed-form rational expressions derived from decoupled Lax integration:
%   s21_exact(lam) = 2 / (1 + i - 2i*lam)^3
%   s11_exact(lam) = 1 + 4 / [(1 + i - 2i*lam)^3 * (1 - i - 2i*lam)^3]
I1_exact  = 2 ./ (1 + 1i - 2i*lambda_vec).^3;
I2_exact  = 2 ./ (1 - 1i - 2i*lambda_vec).^3;
s21_exact = I1_exact;
s11_exact = 1 + I1_exact .* I2_exact;

err_s11 = abs(s11_num - s11_exact);
err_s21 = abs(s21_num - s21_exact);

%% 4. Generate LaTeX Table Code (reflection_table.tex)
tex_filename = 'reflection_table.tex';
fid = fopen(tex_filename, 'w');

tex_lines = {
'\begin{table}[htbp]'
'\centering'
'\small'
'\renewcommand{\arraystretch}{1.15}'
'\caption{Validation of direct scattering coefficients against closed-form expressions for the benchmark profile $q_0(x) = x^2\mathrm{e}^{-(1-\mathrm{i})x}\Theta(x)$, where $\Theta(x)$ denotes the Heaviside step function. Absolute errors are defined as $\mathcal{E}(s) = |s^{\mathrm{num}} - s^{\mathrm{exact}}|$.}'
'\label{tab:reflection_table}'
'\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}ccc}'
'\toprule'
'Spectral parameter $\lambda$ & Absolute error $\mathcal{E}(s_{11})$ & Absolute error $\mathcal{E}(s_{21})$ \\'
'\midrule'
};

for k = 1:num_lambda
    lam = lambda_vec(k);
    if lam >= 0
        lam_str = sprintf('$\\phantom{-}%.2f$', lam);
    else
        lam_str = sprintf('$%.2f$', lam);
    end
    
    e11_str = format_sci_latex(err_s11(k));
    e21_str = format_sci_latex(err_s21(k));
    
    row_str = sprintf('%s & %s & %s \\\\', lam_str, e11_str, e21_str);
    tex_lines{end+1} = row_str;
end

tex_lines{end+1} = '\bottomrule';
tex_lines{end+1} = '\end{tabular*}';
tex_lines{end+1} = '\end{table}';

% Display in command window and export to file
fprintf('\n================== Generated LaTeX Table ==================\n\n');
for i = 1:length(tex_lines)
    fprintf('%s\n', tex_lines{i});
    if fid > 0
        fprintf(fid, '%s\n', tex_lines{i});
    end
end
if fid > 0
    fclose(fid);
    fprintf('\n>> LaTeX table successfully exported to: %s\n', tex_filename);
end
fprintf('===========================================================\n');

%% 5. Display Summary Statistics
fprintf('\n================== Accuracy Statistics ====================\n');
fprintf('s11 Error Summary:\n');
fprintf('  Max Error:  %.4e (at lambda = %+.2f)\n', max(err_s11), lambda_vec(err_s11 == max(err_s11)));
fprintf('  Mean Error: %.4e\n', mean(err_s11));
fprintf('s21 Error Summary:\n');
fprintf('  Max Error:  %.4e (at lambda = %+.2f)\n', max(err_s21), lambda_vec(err_s21 == max(err_s21)));
fprintf('  Mean Error: %.4e\n', mean(err_s21));
fprintf('===========================================================\n');

%% =========================================================================
%% Helper Function: Format numerical scalars into LaTeX scientific notation
%% =========================================================================
function str = format_sci_latex(val)
    if val == 0 || abs(val) < 1e-15
        str = '$0$';
        return;
    end
    s = sprintf('%.4e', val);
    parts = strsplit(s, 'e');
    mantissa = str2double(parts{1});
    exponent = str2double(parts{2});
    str = sprintf('$%.4f \\times 10^{%d}$', mantissa, exponent);
end