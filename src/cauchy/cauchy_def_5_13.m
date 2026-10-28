function C = cauchy_def_5_13(x, n, def_src, def_tgt)
%CAUCHY_DEF_5_13 Evaluate Cauchy operator for ray source and disjoint target.
%   C = CAUCHY_DEF_5_13(X, N, DEF_SRC, DEF_TGT) evaluates the discrete
%   Cauchy transform matrix from an unbounded ray to a disjoint target contour.
%
%   Inputs:
%       x       - Collocation nodes on target curve (length n_Omega)
%       n       - Number of collocation nodes on source ray
%       def_src - Source ray definition structure
%       def_tgt - Disjoint target contour definition structure
%
%   Outputs:
%       C       - Discrete Cauchy transform matrix of size length(x) x n

    m = n + 1;

    if strcmp(def_tgt.type, 'ray')
        if isfield(def_tgt, 'direction') && strcmp(def_tgt.direction, 'backward')
            x_tgt_full = [-1; x(:)];
        else
            x_tgt_full = [x(:); 1];
        end
    else
        x_tgt_full = x(:);
    end

    pts_tgt = curve_mapping(def_tgt, x_tgt_full);
    z_omega = pts_tgt{1};
    w = forward_map(def_src, z_omega);

    z_psi = w - sqrt(w - 1) .* sqrt(1 + w);
    Psi_z = psi_row(z_psi, m);

    x_I = flipud(cos(pi * (0:m-1)' / (m - 1)));
    F = inv(cheby(x_I, m));
    C_base = Psi_z * F;

    if strcmp(def_src.direction, 'forward')
        M_inf = 1;
    else
        M_inf = -1;
    end
    [mu_L, mu_R] = algorithm_5_4(m);

    if M_inf == -1
        C_sliced = C_base(:, 2:m);
        F_gamma  = F(:, 2:m);
        mu_vec   = mu_L;
    else
        C_sliced = C_base(:, 1:n);
        F_gamma  = F(:, 1:n);
        mu_vec   = mu_R;
    end

    sub = ones(size(C_sliced, 1), 1) * (mu_vec * F_gamma);
    C = C_sliced - sub;
end