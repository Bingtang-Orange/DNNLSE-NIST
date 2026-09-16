function T = cheby(x, n)
% CHEBY Evaluates the first n Chebyshev polynomials T_0(x), ..., T_{n-1}(x).
%
% Syntax:
%   T = cheby(x, n)
%
% Inputs:
%   x - Evaluation coordinates (column vector; real or complex)
%   n - Polynomial basis dimension (positive integer)
%
% Outputs:
%   T - Matrix of size length(x) x n, where the k-th column corresponds to T_{k-1}(x)
%
% Description:
%   Evaluates Chebyshev polynomials via the three-term recurrence relation:
%       T_0(x) = 1,  T_1(x) = x,  T_k(x) = 2*x*T_{k-1}(x) - T_{k-2}(x).
%   The algebraic recurrence is numerically stable across the complex plane
%   and avoids trigonometric precision loss near the interval endpoints.
%
% References:
%   L. N. Trefethen, "Spectral Methods in MATLAB", SIAM, Philadelphia (2000).

    x = x(:);
    N = length(x);
    T = zeros(N, n);

    if n >= 1
        T(:, 1) = 1;
    end
    if n >= 2
        T(:, 2) = x;
    end
    for k = 3:n
        T(:, k) = 2 * x .* T(:, k-1) - T(:, k-2);
    end
end