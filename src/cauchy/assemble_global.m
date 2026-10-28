function [Cplus, Cminus] = assemble_global(defs, x_all, n_all)
%ASSEMBLE_GLOBAL Assemble global Cauchy operator matrices C^+ and C^-.
%   [CPLUS, CMINUS] = ASSEMBLE_GLOBAL(DEFS, X_ALL, N_ALL) constructs the
%   discrete boundary Cauchy transform operators on the full contour system.
%
%   Inputs:
%       defs   - Cell array of contour definitions
%       x_all  - Concatenated collocation nodes on [-1, 1] across all contours
%       n_all  - Vector of collocation node counts for each contour
%
%   Outputs:
%       Cplus  - Global discrete boundary Cauchy operator C^+ (size N x N)
%       Cminus - Global discrete boundary Cauchy operator C^- (size N x N)

    num_curves = length(defs);
    N = sum(n_all);
    Cplus = zeros(N, N);

    block_start = cumsum([0, n_all(1:end-1)]) + 1;
    block_end   = cumsum(n_all);

    for i = 1:num_curves
        idx_tgt = block_start(i):block_end(i);
        x_tgt   = x_all(idx_tgt);
        for j = 1:num_curves
            idx_src = block_start(j):block_end(j);
            is_self = (i == j);
            Cplus(idx_tgt, idx_src) = assemble_cauchy_block(x_tgt, n_all(j), defs{j}, defs{i}, is_self);
        end
    end

    Cminus = Cplus;
    for k = 1:num_curves
        idx = block_start(k):block_end(k);
        Cminus(idx, idx) = Cplus(idx, idx) - eye(n_all(k));
    end
end