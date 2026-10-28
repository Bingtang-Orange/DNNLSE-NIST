function C = cauchy_def_5_3(x, n_I, def)
%CAUCHY_DEF_5_3 Evaluate Cauchy matrix from the unit interval to a disjoint contour.
%   C = CAUCHY_DEF_5_3(X, N_I, DEF) computes the discrete Cauchy transform
%   from the reference interval [-1, 1] to a disjoint exterior contour.
%
%   Inputs:
%       x   - Reference collocation nodes mapped from the target curve
%       n_I - Polynomial basis dimension on the source interval
%       def - Geometric structure defining the target contour
%
%   Outputs:
%       C   - Discrete Cauchy transform matrix of size length(x) x n_I

    if strcmp(def.type, 'ray')
        if isfield(def, 'direction') && ~isempty(def.direction) && strcmp(def.direction, 'backward')
            x_full = [-1; x];
        else
            x_full = [x; 1];
        end
    else
        x_full = x;
    end

    pts = curve_mapping(def, x_full);
    z_curve = pts{1};

    z_psi = z_curve - sqrt(z_curve - 1) .* sqrt(1 + z_curve);

    Psi_z = psi_row(z_psi, n_I);

    x_I = flipud(cos(pi*(0:n_I-1)'/(n_I-1)));
    F = inv(cheby(x_I, n_I));

    C = Psi_z * F;
end