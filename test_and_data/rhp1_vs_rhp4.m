% Unit test: Consistency check between RHP I and RHP IV at t = 1.0.

clear; clc;

%% Parameters and potential
q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-15, 15]);
t_val = 1.0;

%% Direct scattering data for RHP I
L_rhp1 = 5.0;
n_rhp1 = 400;
x_I_1 = flipud(cos(pi * (0:n_rhp1-1)' / (n_rhp1 - 1)));
def_1 = create_curve_def('interval', -L_rhp1, L_rhp1);
pts_1 = curve_mapping(def_1, x_I_1);
lambda_1 = pts_1{1};

[rho_1, rho_hat_1] = compute_reflection(lambda_1, q0, L_rhp1);

%% Solver parameters and test nodes
r_box_t1 = 0.25;
d_ray_t1 = 2.20;
n_box    = 24;
n_ray    = 24;

x_test = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];

%% Compute solutions and compare
fprintf('========================================================================================\n');
fprintf('  Consistency Check: RHP I vs RHP IV (t = 1.0)\n');
fprintf('========================================================================================\n');
fprintf('%-6s | %-8s | %-24s | %-24s | %-10s | %-12s\n', ...
        'x', 'lambda0', 'q (RHP I)', 'q (RHP IV)', '|Diff|', 'Zero-Sum Res');
fprintf('----------------------------------------------------------------------------------------\n');

for j = 1:length(x_test)
    xev = x_test(j);
    
    q1 = solve_rhp1(xev, t_val, lambda_1, rho_1, rho_hat_1, L_rhp1, n_rhp1);
    [q4, info4] = solve_rhp4(xev, t_val, q0, 8.0, 8.0, n_box, n_ray, r_box_t1, d_ray_t1);
    
    err_diff = abs(q4 - q1);
    
    fprintf('%+6.2f | %+8.3f | %+.7f%+.7fi | %+.7f%+.7fi | %.3e | %.3e\n', ...
            xev, info4.lambda0, real(q1), imag(q1), real(q4), imag(q4), err_diff, info4.zero_sum_max);
end
fprintf('----------------------------------------------------------------------------------------\n');