function C = cauchy_def_5_1(x, n)
%CAUCHY_DEF_5_1 Evaluate discrete Cauchy operator C^+ on the unit interval.
%   C = CAUCHY_DEF_5_1(X, N) computes the discrete boundary Cauchy matrix
%   C^+ on the reference interval [-1, 1].
%
%   Inputs:
%       x - Collocation nodes on [-1, 1] (including endpoints)
%       n - Number of Chebyshev collocation nodes
%
%   Outputs:
%       C - Discrete Cauchy matrix of size n x n

    [C, ~, ~, ~, ~] = cauchy_kernel_interval(x, n);
end