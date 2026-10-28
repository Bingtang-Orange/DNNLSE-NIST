function T = cheby(x, n)
%CHEBY Evaluate Chebyshev polynomials of the first kind.
%   T = CHEBY(X, N) computes values of Chebyshev polynomials T_0(x), ..., T_{n-1}(x)
%   at points X using the three-term recurrence relation.
%
%   Inputs:
%       x - Evaluation coordinates (column vector)
%       n - Polynomial basis dimension (positive integer)
%
%   Outputs:
%       T - Matrix of size length(x) x n, where column k corresponds to T_{k-1}(x)

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