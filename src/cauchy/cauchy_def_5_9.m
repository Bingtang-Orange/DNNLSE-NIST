function C_Gamma = cauchy_def_5_9(x, n, def)
% CAUCHY_DEF_5_9 Discrete Cauchy operator C^+[Gamma, Gamma] for arcs (M_Gamma(inf) ~= inf).
%
% Syntax:
%   C_Gamma = cauchy_def_5_9(x, n, def)
%
% Inputs:
%   x       - Collocation nodes on the reference unit interval
%   n       - Number of Chebyshev collocation nodes
%   def     - Bounded arc definition struct containing finite field M_inf
%
% Outputs:
%   C_Gamma - Modified discrete Cauchy operator matrix on Gamma
%
% Description:
%   For circular arcs where M_Gamma(inf) is finite, the standard operator
%   cauchy_def_5_6 is regularized by subtracting a rank-one projection matrix
%   evaluated at the conformal image of infinity z = T_+^{-1}(M_inf).
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.9.

    % 1. Evaluate base diagonal endpoint regularized operator (Definition 5.6)
    C_Gamma = cauchy_def_5_6(x, n, def);

    % 2. Validate finite conformal image of infinity
    if ~isfield(def, 'M_inf') || isempty(def.M_inf) || isinf(def.M_inf)
        error('cauchy_def_5_9: def.M_inf must be a finite scalar. Verify curve definition.');
    end
    M_inf = def.M_inf;

    % 3. Evaluate Joukowsky conformal inverse z = T_+^{-1}(M_inf)
    z = M_inf - sqrt(M_inf - 1) .* sqrt(M_inf + 1);

    % 4. Construct basis evaluation row vector Psi_z
    Psi_z = psi_row(z, n);

    % 5. Chebyshev basis transformation matrix F
    F = inv(cheby(x, n));

    % 6. Subtract rank-one correction matrix: ones(n, 1) * (Psi_z * F)
    sub = ones(n, 1) * (Psi_z * F);
    C_Gamma = C_Gamma - sub;
end