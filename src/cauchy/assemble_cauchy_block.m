function C = assemble_cauchy_block(x_tgt, n_src, def_src, def_tgt, is_self)
%ASSEMBLE_CAUCHY_BLOCK Construct discrete Cauchy transform interaction block.
%   C = ASSEMBLE_CAUCHY_BLOCK(X_TGT, N_SRC, DEF_SRC, DEF_TGT)
%   C = ASSEMBLE_CAUCHY_BLOCK(X_TGT, N_SRC, DEF_SRC, DEF_TGT, IS_SELF)
%   evaluates the discrete Cauchy transform block matrix between source
%   and target contours based on their geometric types and connectivity.
%
%   Inputs:
%       x_tgt   - Collocation nodes on [-1, 1] for the target contour
%       n_src   - Polynomial basis dimension on the source contour
%       def_src - Geometric structure defining the source contour
%       def_tgt - Geometric structure defining the target contour
%       is_self - (Optional) Boolean flag indicating self-interaction
%
%   Outputs:
%       C       - Discrete Cauchy transform block matrix of size length(x_tgt) x n_src

    if nargin < 5 || isempty(is_self)
        is_self = check_identical_curve(def_src, def_tgt);
    end

    is_identity_src = strcmp(def_src.type, 'interval') && ...
                      abs(def_src.a + 1) < 1e-12 && abs(def_src.b - 1) < 1e-12;

    if is_self
        is_connected = true;
    else
        is_connected = check_connected(def_src, def_tgt);
    end

    if is_identity_src
        if is_self
            C = cauchy_def_5_1(x_tgt, n_src);
        elseif is_connected
            C = cauchy_def_5_5(x_tgt, n_src, def_tgt);
        else
            C = cauchy_def_5_3(x_tgt, n_src, def_tgt);
        end
    elseif strcmp(def_src.type, 'interval')
        if is_self
            C = cauchy_def_5_6(x_tgt, n_src, def_src);
        elseif is_connected
            C = cauchy_def_5_8(x_tgt, n_src, def_src, def_tgt);
        else
            C = cauchy_def_5_7(x_tgt, n_src, def_src, def_tgt);
        end
    elseif strcmp(def_src.type, 'arc')
        if is_self
            C = cauchy_def_5_9(x_tgt, n_src, def_src);
        elseif is_connected
            C = cauchy_def_5_11(x_tgt, n_src, def_src, def_tgt);
        else
            C = cauchy_def_5_10(x_tgt, n_src, def_src, def_tgt);
        end
    elseif strcmp(def_src.type, 'ray')
        if is_self
            C = cauchy_def_5_12(x_tgt, n_src, def_src);
        elseif is_connected
            C = cauchy_def_5_14(x_tgt, n_src, def_src, def_tgt);
        else
            C = cauchy_def_5_13(x_tgt, n_src, def_src, def_tgt);
        end
    else
        error('assemble_cauchy_block: Unsupported source contour type: %s', def_src.type);
    end
end

function flag = check_connected(def_src, def_tgt)
    tol = 1e-11;
    src_ends = get_ends(def_src);
    tgt_ends = get_ends(def_tgt);
    for i = 1:length(src_ends)
        for j = 1:length(tgt_ends)
            if isfinite(src_ends(i)) && isfinite(tgt_ends(j))
                if abs(src_ends(i) - tgt_ends(j)) < tol
                    flag = true;
                    return;
                end
            end
        end
    end
    flag = false;
end

function flag = check_identical_curve(d1, d2)
    tol = 1e-12;
    if ~strcmp(d1.type, d2.type)
        flag = false;
        return;
    end
    e1 = get_ends(d1);
    e2 = get_ends(d2);
    flag = all(abs(e1 - e2) < tol);
end

function ends = get_ends(def)
    switch def.type
        case 'interval'
            ends = [def.a, def.b];
        case 'ray'
            dir = 'forward';
            if isfield(def, 'direction') && ~isempty(def.direction)
                dir = def.direction;
            end
            if strcmp(dir, 'forward')
                ends = [def.a, Inf];
            else
                ends = [Inf, def.a];
            end
        case 'arc'
            ends = [def.a + def.r*exp(1i*def.theta1), def.a + def.r*exp(1i*def.theta2)];
    end
end