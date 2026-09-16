% RUN_STEP_VALIDATION_TABLE
% Benchmarks Chebfun direct scattering against exact/quadrature reference data.
% Validates s_11 and s_21 for step-like initial potential q(x) = x^2 * exp(-x + ix) (x >= 0).
% Reproduces Section 3.2.1 verification results.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

% ---- Configuration ----
L = 40.0;
lambda_vec = -2.0:0.5:2.0;
num_lambda = length(lambda_vec);

% Step-like initial potential
q_step = @(x) (x >= 0) .* (x.^2 .* exp(-x) .* exp(1i*x));
q = chebfun(q_step, [-L, L], 'splitting', 'on');

% Prepare potentials struct for fast solver execution
x_left  = chebfun('x', [-L, 0]);
x_right = chebfun('x', [0, L]);
pot.q_left        = q(x_left);
pot.q_neg_c_left  = conj(q(-x_left));
pot.q_right       = q(x_right);
pot.q_neg_c_right = conj(q(-x_right));
pot.L             = L;

%% 1. Compute Numerical Solutions via Chebfun
s11_num = complex(zeros(1, num_lambda));
s21_num = complex(zeros(1, num_lambda));

fprintf('Computing Chebfun ODE solutions on spectral grid...\n');
for k = 1:num_lambda
    S = compute_S(lambda_vec(k), pot);
    s11_num(k) = S(1,1);
    s21_num(k) = S(2,1);
end

%% 2. Compute Analytical Reference Values via Exact Formulas & Adaptive Quadrature
s11_ref = complex(zeros(1, num_lambda));
s21_ref = complex(zeros(1, num_lambda));

fprintf('Evaluating reference integrals...\n');
for k = 1:num_lambda
    lam = lambda_vec(k);
    % Exact closed-form formula (Eq. (26) in Section 3.2.1)
    s11_ref(k) = 1 + 4 / (((1 + 1i - 2i*lam)^3) * ((1 - 1i - 2i*lam)^3));
    s21_ref(k) = 2 / ((1 + 1i - 2i*lam)^3);
end

%% 3. Print Validation Table
err_s11 = abs(s11_num - s11_ref);
err_s21 = abs(s21_num - s21_ref);

fprintf('\n=======================================================================================================\n');
fprintf(' Validation Table for Step-like Initial Data (q(x) = x^2 exp(-x+ix), x >= 0)\n');
fprintf('=======================================================================================================\n');
fprintf('%-6s | %-22s %-22s %-12s | %-12s\n', ...
        'lambda', 's11 (Chebfun)', 's11 (Exact)', 'Err(s11)', 'Err(s21)');
fprintf('-------------------------------------------------------------------------------------------------------\n');

for k = 1:num_lambda
    s11_str = sprintf('%.4e%+.4ei', real(s11_num(k)), imag(s11_num(k)));
    ref_str = sprintf('%.4e%+.4ei', real(s11_ref(k)), imag(s11_ref(k)));
    fprintf('%-6.2f | %-22s %-22s %-12.3e | %-12.3e\n', ...
            lambda_vec(k), s11_str, ref_str, err_s11(k), err_s21(k));
end
fprintf('-------------------------------------------------------------------------------------------------------\n');

% Summary statistics
fprintf('\nSummary Statistics:\n');
fprintf('  Max Error in s11: %.3e (Mean: %.3e)\n', max(err_s11), mean(err_s11));
fprintf('  Max Error in s21: %.3e (Mean: %.3e)\n', max(err_s21), mean(err_s21));
fprintf('=======================================================================================================\n');