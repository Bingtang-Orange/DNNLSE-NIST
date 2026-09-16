function C = cauchy_def_5_3(x, n_I, def)
% CAUCHY_DEF_5_3 Computes Cauchy matrix C[I, Gamma] from unit interval to disjoint curve.
%
% Syntax:
%   C = cauchy_def_5_3(x, n_I, def)
%
% Inputs:
%   x    - Reference collocation nodes mapped from the target curve Gamma (length n_Gamma)
%   n_I  - Polynomial basis dimension on the source interval I = [-1, 1]
%   def  - Geometric struct defining the target curve Gamma (via create_curve_def)
%
% Outputs:
%   C    - Discrete Cauchy transform matrix C[I, Gamma] of size n_Gamma x n_I
%
% Description:
%   Evaluates the Cauchy integral on an exterior disjoint contour Gamma (where
%   Gamma \cap I = \emptyset) using the inverse Joukowsky-type conformal map
%   T_+^{-1}(z) = z - sqrt(z - 1)*sqrt(z + 1), mapping the exterior of [-1, 1]
%   into the interior of the unit disk.
%
% References:
%   S. Olver, Numer. Math. (2012), Definition 5.3.

    % 1. Restore point at infinity for unbounded rays prior to domain slicing
    if strcmp(def.type, 'ray')
        if isfield(def, 'direction') && ~isempty(def.direction) && strcmp(def.direction, 'backward')
            % Backward ray: point at infinity mapped to x = -1
            x_full = [-1; x];
        else
            % Forward ray: point at infinity mapped to x = +1
            x_full = [x; 1];
        end
    else
        % Bounded contours (intervals, arcs): no truncation slicing required
        x_full = x;
    end

    % 2. Map mapped nodes to target contour coordinates in the complex plane
    pts = curve_mapping(def, x_full);
    z_curve = pts{1};

    % 3. Evaluate conformal mapping T_+^{-1}(z)
    z_psi = z_curve - sqrt(z_curve - 1) .* sqrt(1 + z_curve);

    % 4. Construct basis evaluations Psi_{n_I}^z of size length(x) x n_I
    Psi_z = psi_row(z_psi, n_I);

    % 5. Compute Chebyshev transformation matrix F on the reference interval
    x_I = flipud(cos(pi*(0:n_I-1)'/(n_I-1)));
    F = inv(cheby(x_I, n_I));

    % 6. Assemble the discrete operator matrix C[I, Gamma]
    C = Psi_z * F;
end