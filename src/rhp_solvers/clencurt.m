function w = clencurt(n)
% CLENCURT Computes Clenshaw-Curtis quadrature weights on Chebyshev-Lobatto nodes.
%
% Syntax:
%   w = clencurt(n)
%
% Inputs:
%   n - Number of quadrature nodes in [-1, 1] including endpoints (n >= 1)
%
% Outputs:
%   w - Row vector of quadrature weights of length n, satisfying:
%       \int_{-1}^{1} f(x) dx \approx w * f(x)
%
% Description:
%   Evaluates exact Clenshaw-Curtis integration weights for Chebyshev-Lobatto 
%   collocation points x_j = cos(j*pi/(n-1)). Used for high-precision numerical
%   quadrature in scalar delta(lambda) evaluations and asymptotic potential reconstruction.
%
% References:
%   C. W. Clenshaw and A. R. Curtis, Numer. Math. 2 (1960), pp. 197-205.
%   L. N. Trefethen, "Spectral Methods in MATLAB", SIAM, Philadelphia (2000), Chapter 12.

    N = n - 1; % Polynomial degree
    if N == 0
        w = 2;
        return;
    end

    theta = pi * (0:N)' / N;
    w = zeros(1, N+1);
    ii = 2:N; % Interior node indices
    v = ones(N-1, 1);

    if mod(N, 2) == 0
        w(1)   = 1 / (N^2 - 1);
        w(N+1) = w(1);
        for k = 1 : (N/2 - 1)
            v = v - 2 * cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
        v = v - cos(N*theta(ii)) / (N^2 - 1);
    else
        w(1)   = 1 / N^2;
        w(N+1) = w(1);
        for k = 1 : ((N-1)/2)
            v = v - 2 * cos(2*k*theta(ii)) / (4*k^2 - 1);
        end
    end

    w(ii) = 2 * v / N;
end