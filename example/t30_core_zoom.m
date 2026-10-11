% Comparison of numerical solution and asymptotic prediction in core region at t = 30.

clear; clc; close all;

%% Load data and setup
load('t30_rhp4_data.mat');

idx_core = (x_grid >= -50.0) & (x_grid <= 50.0);
x_core   = x_grid(idx_core);
q_core   = q_vals(idx_core);

q0_handle = @(x) 0.9 * exp(-x.^2 + x);
q0 = chebfun(q0_handle, [-20, 20]);

%% Precompute reflection spectrum
s_dense = linspace(-6.0, 6.0, 601).';
ds = s_dense(2) - s_dense(1);
r_dense  = complex(zeros(size(s_dense)));
rh_dense = complex(zeros(size(s_dense)));

for k = 1:length(s_dense)
    [r_dense(k), rh_dense(k)] = compute_reflection(s_dense(k), q0, 8.0);
end
ln_det   = log(1.0 - r_dense .* rh_dense);
d_ln_det = gradient(ln_det, ds);

%% Evaluate asymptotic discrete samples
dx_asymp  = 1.0; 
x_asymp   = (-50.0:dx_asymp:50.0).';
num_asymp = length(x_asymp);
q_asymp   = complex(zeros(num_asymp, 1));

for j = 1:num_asymp
    xj = x_asymp(j);
    lam0 = xj / (4.0 * t_val);
    
    [r0, rh0] = compute_reflection(lam0, q0, 8.0);
    det0 = 1.0 - r0 * rh0;
    nu0  = -(1.0 / (2.0 * pi)) * log(det0);
    nu_R = real(nu0);
    nu_I = imag(nu0);
    
    idx_tail = find(s_dense <= lam0);
    if length(idx_tail) >= 3
        s_sub = s_dense(idx_tail);
        u = lam0 - s_sub;
        u(end) = max(u(end), 1e-6);
        kernel = log(u) .* d_ln_det(idx_tail);
        chi_val = -(1.0 / (2.0 * pi * 1i)) * trapz(s_sub, kernel);
    else
        chi_val = 0.0;
    end
    
    scale_fac = exp(-3.0 * 1i * nu0 * log(2.0)); 
    gamma_val = cgamma_lanczos(-1i * nu0);
    
    numerator = sqrt(pi) * exp(-pi*nu0/2.0 + 1i*pi/4.0 + 2.0*chi_val) * scale_fac;
    p_fac     = numerator / (rh0 * gamma_val);
    
    carrier_phase = exp(-1i * (xj^2) / (4.0 * t_val));
    log_factor    = exp(1i * nu_R * log(t_val)); 
    decay_factor  = t_val^(-0.50 - nu_I);
    
    q_asymp(j) = - decay_factor * p_fac * carrier_phase * log_factor;
end

%% Plot comparison in core region
figure('Name', 'Core Region Zoom', 'Color', 'w', 'Position', [150, 150, 1050, 650]);

subplot(2, 1, 1);
plot(x_core, real(q_core), 'r-', 'LineWidth', 1.1, 'DisplayName', 'Numerical Re q'); hold on;
plot(x_asymp, real(q_asymp), 'ko', 'MarkerSize', 4.0, 'LineWidth', 0.8, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic Re q');
grid on; box on;
xlim([-50, 50]);
ylim([-0.15, 0.15]);
ylabel('$\mathrm{Re}\,q(x, 30)$', 'Interpreter', 'latex', 'FontSize', 12);
title('Core Region Zoom: $x \in [-50, 50]$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');

subplot(2, 1, 2);
plot(x_core, imag(q_core), 'b-', 'LineWidth', 1.1, 'DisplayName', 'Numerical Im q'); hold on;
plot(x_asymp, imag(q_asymp), 'ko', 'MarkerSize', 4.0, 'LineWidth', 0.8, ...
     'MarkerFaceColor', 'none', 'DisplayName', 'Asymptotic Im q');
grid on; box on;
xlim([-50, 50]);
ylim([-0.15, 0.15]);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\mathrm{Im}\,q(x, 30)$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');

%% Helper functions
function g = cgamma_lanczos(z)
    p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, ...
         771.32342877765313, -176.61502916214059, 12.507343278686905, ...
         -0.138571095836524, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if real(z) < 0.5
        g = pi ./ (sin(pi .* z) .* cgamma_lanczos(1.0 - z));
    else
        z_adj = z - 1.0;
        x = p(1);
        for k = 2:length(p)
            x = x + p(k) ./ (z_adj + (k - 1));
        end
        t = z_adj + 7.5;
        g = sqrt(2.0 * pi) .* (t.^(z_adj + 0.5)) .* exp(-t) .* x;
    end
end