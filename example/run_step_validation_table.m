clear; clc;
warning('off', 'chebfun:constructor');
warning('off', 'chebfun:compose');

L = 40;
lambda_vec = -2.0:0.5:2.0;
num_lambda = length(lambda_vec);

f = @(s) s.^2 .* exp(-s) .* exp(1i*s);

pot = struct();
pot.L             = L;
pot.q_left        = chebfun(0, [-L, 0]);
pot.q_neg_c_left  = chebfun(@(x) conj(f(-x)), [-L, 0], 'splitting', 'on');
pot.q_right       = chebfun(@(x) f(x), [0, L], 'splitting', 'on');
pot.q_neg_c_right = chebfun(0, [0, L]);

s11_num = zeros(1, num_lambda);
s21_num = zeros(1, num_lambda);

fprintf('Computing numerical scattering data via Chebyshev collocation...\n');
for k = 1:num_lambda
    lam = lambda_vec(k);
    S = compute_S(lam, pot);
    s11_num(k) = S(1,1);
    s21_num(k) = S(2,1);
end

I1_exact  = 2 ./ (1 + 1i - 2i*lambda_vec).^3;
I2_exact  = 2 ./ (1 - 1i - 2i*lambda_vec).^3;
s21_exact = I1_exact;
s11_exact = 1 + I1_exact .* I2_exact;

err_s11 = abs(s11_num - s11_exact);
err_s21 = abs(s21_num - s21_exact);

fprintf('\n===================================================================\n');
fprintf('  Validation of Direct Scattering Coefficients\n');
fprintf('===================================================================\n');
fprintf('  %-10s | %-24s | %-24s\n', 'lambda', 'Error of s11', 'Error of s21');
fprintf('-------------------------------------------------------------------\n');
for k = 1:num_lambda
    fprintf('  %+6.2f     | %18.4e       | %18.4e\n', ...
            lambda_vec(k), err_s11(k), err_s21(k));
end
fprintf('===================================================================\n\n');

fprintf('================== Accuracy Statistics ====================\n');
fprintf('s11 Error Summary:\n');
fprintf('  Max Error:  %.4e (at lambda = %+.2f)\n', max(err_s11), lambda_vec(err_s11 == max(err_s11)));
fprintf('  Mean Error: %.4e\n', mean(err_s11));
fprintf('s21 Error Summary:\n');
fprintf('  Max Error:  %.4e (at lambda = %+.2f)\n', max(err_s21), lambda_vec(err_s21 == max(err_s21)));
fprintf('  Mean Error: %.4e\n', mean(err_s21));
fprintf('===========================================================\n');