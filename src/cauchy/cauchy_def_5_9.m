function C_Gamma = cauchy_def_5_9(x, n, def)
%CAUCHY_DEF_5_9 Discrete Cauchy operator C^+ on circular arcs.
%   C_GAMMA = CAUCHY_DEF_5_9(X, N, DEF) computes the discrete Cauchy
%   operator matrix on circular arcs.
%
%   Inputs:
%       x       - Collocation nodes on the reference unit interval
%       n       - Number of Chebyshev collocation nodes
%       def     - Bounded arc definition structure containing scalar M_inf
%
%   Outputs:
%       C_Gamma - Modified discrete Cauchy operator matrix on the arc

    C_Gamma = cauchy_def_5_6(x, n, def);

    if ~isfield(def, 'M_inf') || isempty(def.M_inf) || isinf(def.M_inf)
        error('cauchy_def_5_9: def.M_inf must be a finite scalar. Verify curve definition.');
    end
    M_inf = def.M_inf;

    z = M_inf - sqrt(M_inf - 1) .* sqrt(M_inf + 1);
    Psi_z = psi_row(z, n);

    F = inv(cheby(x, n));
    sub = ones(n, 1) * (Psi_z * F);
    C_Gamma = C_Gamma - sub;
end