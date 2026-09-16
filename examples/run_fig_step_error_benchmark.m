% RUN_FIG_STEP_ERROR_BENCHMARK
% Benchmarks Chebfun direct scattering against exact closed-form solutions
% across four step-like initial potentials with varying boundary regularities.
% Reproduces Figure 4 in the Response to Reviewers and Section 3.2.1 results.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

%% 1. Simulation Parameters and Benchmark Potentials
L = 40.0;                    % Truncation length ensuring machine precision
lambda_vec = -2.0:0.5:2.0;   % Spectral parameter sampling grid
num_lambda = length(lambda_vec);
ICs = struct();

% Profile 1: f(s) = s^2 e^{-s} (C^1 at origin)
ICs(1).name   = 's^2 e^{-s}';
ICs(1).label  = '$s^2 \mathrm{e}^{-s}$';
ICs(1).f      = @(s) s.^2 .* exp(-s);
ICs(1).I1_expr= @(lam) 2 ./ (1 - 2i*lam).^3;
ICs(1).I2_expr= @(lam) 2 ./ (1 - 2i*lam).^3;

% Profile 2: f(s) = s^2 e^{-s} e^{is} (Complex-valued)
ICs(2).name   = 's^2 e^{-s} e^{is}';
ICs(2).label  = '$s^2 \mathrm{e}^{-s}\mathrm{e}^{\mathrm{i}s}$';
ICs(2).f      = @(s) s.^2 .* exp(-s) .* exp(1i*s);
ICs(2).I1_expr= @(lam) 2 ./ (1 + 1i - 2i*lam).^3;
ICs(2).I2_expr= @(lam) 2 ./ (1 - 1i - 2i*lam).^3;

% Profile 3: f(s) = s e^{-s} (C^0 at origin)
ICs(3).name   = 's e^{-s}';
ICs(3).label  = '$s \mathrm{e}^{-s}$';
ICs(3).f      = @(s) s .* exp(-s);
ICs(3).I1_expr= @(lam) 1 ./ (1 - 2i*lam).^2;
ICs(3).I2_expr= @(lam) 1 ./ (1 - 2i*lam).^2;

% Profile 4: f(s) = s^3 e^{-s} (C^2 at origin)
ICs(4).name   = 's^3 e^{-s}';
ICs(4).label  = '$s^3 \mathrm{e}^{-s}$';
ICs(4).f      = @(s) s.^3 .* exp(-s);
ICs(4).I1_expr= @(lam) 6 ./ (1 - 2i*lam).^4;
ICs(4).I2_expr= @(lam) 6 ./ (1 - 2i*lam).^4;

num_ICs = length(ICs);
err_s11 = zeros(num_ICs, num_lambda);
err_s21 = zeros(num_ICs, num_lambda);

fprintf('=========================================================================\n');
fprintf(' Pointwise Verification against Closed-Form Exact Scattering Data\n');
fprintf('=========================================================================\n');

%% 2. Direct Scattering Computation Loop
for idx = 1:num_ICs
    f = ICs(idx).f;
    fprintf('--> Computing for IC: %-20s ... ', ICs(idx).name);
    
    % Construct Chebfun step potential: q(x) = f(x) for x >= 0, 0 for x < 0
    q = chebfun(@(x) (x >= 0) .* f(x), [-L, L], 'splitting', 'on');
    
    % Pre-sample potentials to construct solver struct
    x_left  = chebfun('x', [-L, 0]);
    x_right = chebfun('x', [0, L]);
    pot.q_left        = q(x_left);
    pot.q_neg_c_left  = conj(q(-x_left));
    pot.q_right       = q(x_right);
    pot.q_neg_c_right = conj(q(-x_right));
    pot.L             = L;
    
    s11_arr = complex(zeros(1, num_lambda));
    s21_arr = complex(zeros(1, num_lambda));
    
    % Solve ODE system across the spectral parameter grid
    for k = 1:num_lambda
        S = compute_S(lambda_vec(k), pot);
        s11_arr(k) = S(1,1);
        s21_arr(k) = S(2,1);
    end
    
    % Evaluate analytical closed-form reference values
    I1_ref = ICs(idx).I1_expr(lambda_vec);
    I2_ref = ICs(idx).I2_expr(lambda_vec);
    a_ref  = 1 + I1_ref .* I2_ref;
    
    % Compute absolute discrepancies
    err_s11(idx, :) = abs(s11_arr - a_ref);
    err_s21(idx, :) = abs(s21_arr - I1_ref);
    
    fprintf('Done (Max Err s11: %.2e, s21: %.2e)\n', max(err_s11(idx,:)), max(err_s21(idx,:)));
end

%% 3. Paper Figure: Scatter Error Distributions (Figure 4)
figure('Color', 'w', 'Position', [150, 100, 750, 600]);
colors = lines(num_ICs);

% Panel 1: s11 Error vs. lambda
subplot(2, 1, 1);
hold on;
for idx = 1:num_ICs
    scatter(lambda_vec, err_s11(idx, :), 50, colors(idx,:), 'filled', ...
            'DisplayName', ICs(idx).label);
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 11, 'LineWidth', 1.1);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$s_{11}$ Error', 'Interpreter', 'latex', 'FontSize', 13);
title('$s_{11}$ Pointwise Absolute Error Across Benchmark Profiles', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
ylim([1e-12, 1e-8]);
grid on; box on;

% Panel 2: s21 Error vs. lambda
subplot(2, 1, 2);
hold on;
for idx = 1:num_ICs
    scatter(lambda_vec, err_s21(idx, :), 50, colors(idx,:), 'filled', ...
            'DisplayName', ICs(idx).label);
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 11, 'LineWidth', 1.1);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$s_{21}$ Error', 'Interpreter', 'latex', 'FontSize', 13);
title('$s_{21}$ Pointwise Absolute Error Across Benchmark Profiles', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
ylim([1e-12, 1e-8]);
grid on; box on;

%% 4. Command-Line Verification Summary
fprintf('\n================== Summary of Maximum Benchmark Errors ==================\n');
fprintf('%-22s | %-16s | %-16s\n', 'Initial Potential Profile', 'Max Error (s11)', 'Max Error (s21)');
fprintf('-------------------------------------------------------------------------\n');
for idx = 1:num_ICs
    fprintf('%-25s | %-16.4e | %-16.4e\n', ICs(idx).name, max(err_s11(idx,:)), max(err_s21(idx,:)));
end
fprintf('=========================================================================\n');