function [q_asymp, mod_asymp] = eval_dnnlse_asymptotics(x_vec, t_val, q0)
%EVAL_DNNLSE_ASYMPTOTICS Compute asymptotic solutions of the defocusing nonlocal NLSE.
%   [Q_ASYMP, MOD_ASYMP] = EVAL_DNNLSE_ASYMPTOTICS(X_VEC, T_VAL, Q0) computes
%   the complex asymptotic potential Q_ASYMP and the modulus envelope MOD_ASYMP
%   at spatial coordinates X_VEC and time T_VAL for the initial condition Q0.
%
%   Inputs:
%       X_VEC     - Vector of spatial coordinates
%       T_VAL     - Evaluation time parameter (scalar, T_VAL > 0)
%       Q0        - Chebfun object of the initial condition
%
%   Outputs:
%       Q_ASYMP   - Vector of complex asymptotic values evaluated at (X_VEC, T_VAL)
%       MOD_ASYMP - Vector of absolute values abs(Q_ASYMP)

    num_pts = length(x_vec);
    q_asymp   = complex(zeros(num_pts, 1));
    mod_asymp = zeros(num_pts, 1);
    
    s_dense = linspace(-6.0, 6.0, 601).';
    ds = s_dense(2) - s_dense(1);
    r_dense  = complex(zeros(size(s_dense)));
    rh_dense = complex(zeros(size(s_dense)));
    for k = 1:length(s_dense)
        [r_dense(k), rh_dense(k)] = compute_reflection(s_dense(k), q0, 8.0);
    end
    ln_det = log(1.0 - r_dense .* rh_dense);
    d_ln_det = gradient(ln_det, ds);

    for j = 1:num_pts
        xj = x_vec(j);
        lam0 = xj / (4.0 * t_val);
        
        [r0, rh0] = compute_reflection(lam0, q0, 8.0);
        det0 = 1.0 - r0 * rh0;
        nu0  = -(1.0 / (2.0 * pi)) * log(det0);
        nu_R = real(nu0);
        nu_I = imag(nu0);
        
        idx_tail = find(s_dense >= lam0 + 1e-5);
        if length(idx_tail) >= 3
            s_sub = s_dense(idx_tail);
            kernel = log(s_sub - lam0) .* d_ln_det(idx_tail);
            chi_val = -(1.0 / (2.0 * pi * 1i)) * trapz(s_sub, kernel);
        else
            chi_val = 0.0;
        end
        
        scale_fac = exp(3.0 * 1i * nu0 * log(2.0)); 
        gamma_val = cgamma_lanczos(1i * nu0);
        numerator = sqrt(pi) * exp(-pi*nu0/2.0 - 1i*pi/4.0 + 2.0*chi_val) * scale_fac;
        p_val     = numerator / (rh0 * gamma_val);
        
        carrier_phase = exp(-1i * (xj^2) / (4.0 * t_val));
        log_phase     = exp(1i * nu_R * log(t_val));
        decay_factor  = t_val^(-0.50 - nu_I);
        
        q_val = decay_factor * p_val * carrier_phase * log_phase;
        q_asymp(j)   = q_val;
        mod_asymp(j) = abs(q_val);
    end
end

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