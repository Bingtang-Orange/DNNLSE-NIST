function Psi_z = psi_row(z, n)
% PSI_ROW Evaluates the exterior Cauchy transform basis matrix Psi_n^z.
%
% Syntax:
%   Psi_z = psi_row(z, n)
%
% Inputs:
%   z     - Vector of evaluation points mapped by T_+^{-1} satisfying |z| < 1
%   n     - Number of Chebyshev basis polynomials
%
% Outputs:
%   Psi_z - Matrix of size length(z) x n representing Cauchy basis evaluations
%
% Description:
%   Implements Algorithm 5.2 of Olver (2012). To prevent catastrophic cancellation,
%   a hybrid evaluation scheme is adopted:
%   (i)  For |z| < 0.8: dynamic adaptive Taylor series summation to machine precision;
%   (ii) For |z| >= 0.8: closed-form algebraic evaluation via leading logarithm psi_0.
%   Two-sided finite-difference recurrences generate the intermediate basis sequence.
%
% References:
%   S. Olver, Numer. Math. (2012), Algorithm 5.2.

    z = z(:);
    m = length(z);
    Psi_z = zeros(m, n);

    % Base evaluation psi_0(z)
    psi_0 = (2 / (1i * pi)) * atanh(z);

    if n == 1
        Psi_z(:, 1) = -psi_0;
        return;
    end

    % 1. Evaluate seed term psi_{1-n}(z)
    J = floor((n - 1) / 2);
    p = 2 * (J + 1) - n;

    psi_seed = zeros(m, 1);
    abs_z = abs(z);

    mask_small = (abs_z < 0.8);
    mask_large = ~mask_small;

    % (A) Small-modulus set (|z| < 0.8): adaptive series expansion
    if any(mask_small)
        z_s = z(mask_small);
        w_s = z_s.^2;
        term_s = (2 / (1i * pi)) * (z_s.^p);

        sum_s  = zeros(size(z_s));
        curr_w = ones(size(z_s));
        m_idx  = 0;
        max_iter = 150;
        tol = eps(1.0);

        while true
            denom = 2 * (J + 1 + m_idx) - 1;
            step  = curr_w / denom;
            sum_s = sum_s + step;
            curr_w = curr_w .* w_s;
            m_idx  = m_idx + 1;

            if max(abs(step(:))) < tol || m_idx >= max_iter
                break;
            end
        end
        psi_seed(mask_small) = sum_s .* term_s;
    end

    % (B) Large-modulus set (|z| >= 0.8): closed-form subtraction
    if any(mask_large)
        z_l = z(mask_large);
        mu_l = zeros(size(z_l));
        for j = 1:J
            mu_l = mu_l + (z_l.^(2*j - 1)) / (2*j - 1);
        end
        psi_seed(mask_large) = (z_l.^(1 - n)) .* (psi_0(mask_large) - (2 / (1i * pi)) * mu_l);
    end

    % 2. Full recurrence across indices (Algorithm 5.2)
    psi = zeros(m, 2*n - 1);
    psi(:, n) = psi_0;          % Index k = 0
    psi(:, 1) = psi_seed;       % Index k = 1 - n

    % Negative recurrence: k from (2 - n) down to -1
    for k = (2 - n) : -1
        if mod(k, 2) ~= 0
            psi(:, k + n) = z .* psi(:, k + n - 1) - 2 / (1i * pi * k);
        else
            psi(:, k + n) = z .* psi(:, k + n - 1);
        end
    end

    % Positive recurrence: k from 1 up to (n - 1)
    for k = 1 : (n - 1)
        if mod(k, 2) ~= 0
            psi(:, k + n) = z .* psi(:, k + n - 1) - 2 / (1i * pi * k);
        else
            psi(:, k + n) = z .* psi(:, k + n - 1);
        end
    end

    % 3. Assemble discrete Cauchy basis matrix
    Psi_z(:, 1) = -psi(:, n);
    for i = 2:n
        Psi_z(:, i) = -0.5 * (psi(:, n - i + 1) + psi(:, n + i - 1));
    end
end