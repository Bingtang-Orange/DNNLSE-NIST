clear; clc; close all;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

L = 40;
lambda_vec = -2.0:0.1:2.0;
num_lambda = length(lambda_vec);

ICs = struct();

ICs(1).name    = '$f_1(x)$';
ICs(1).f       = @(s) s .* exp(-s);
ICs(1).I1_expr = @(lam) 1 ./ (1 - 2i*lam).^2;
ICs(1).I2_expr = @(lam) 1 ./ (1 - 2i*lam).^2;

ICs(2).name    = '$f_2(x)$';
ICs(2).f       = @(s) s.^2 .* exp(-s);
ICs(2).I1_expr = @(lam) 2 ./ (1 - 2i*lam).^3;
ICs(2).I2_expr = @(lam) 2 ./ (1 - 2i*lam).^3;

ICs(3).name    = '$f_3(x)$';
ICs(3).f       = @(s) s.^3 .* exp(-s);
ICs(3).I1_expr = @(lam) 6 ./ (1 - 2i*lam).^4;
ICs(3).I2_expr = @(lam) 6 ./ (1 - 2i*lam).^4;

ICs(4).name    = '$f_4(x)$';
ICs(4).f       = @(s) s.^2 .* exp(-s) .* exp(1i*s);
ICs(4).I1_expr = @(lam) 2 ./ (1 + 1i - 2i*lam).^3;
ICs(4).I2_expr = @(lam) 2 ./ (1 - 1i - 2i*lam).^3;

num_ICs = length(ICs);
err_s11 = zeros(num_ICs, num_lambda);
err_s21 = zeros(num_ICs, num_lambda);

for idx = 1:num_ICs
    f = ICs(idx).f;
    
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
    
    I1_ref = ICs(idx).I1_expr(lambda_vec);
    I2_ref = ICs(idx).I2_expr(lambda_vec);
    s11_ref = 1 + I1_ref .* I2_ref;
    s21_ref = I1_ref;
    
    err_s11(idx, :) = abs(s11_arr - s11_ref);
    err_s21(idx, :) = abs(s21_arr - s21_ref);
end

fprintf('\n==================== Error Summary Table ========================\n');
fprintf('%-15s %-18s %-18s\n', 'Benchmark', 'max |Err(s11)|', 'max |Err(s21)|');
fprintf('-----------------------------------------------------------------\n');
for idx = 1:num_ICs
    fprintf('%-15s %-18.4e %-18.4e\n', ICs(idx).name, max(err_s11(idx,:)), max(err_s21(idx,:)));
end
fprintf('=================================================================\n\n');

figure('Color', 'w', 'Position', [100, 100, 800, 600]);
colors  = {'#0072BD', '#D95319', '#2CA02C', '#7E2F8E'};
markers = {'o', 's', '^', 'd'};

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