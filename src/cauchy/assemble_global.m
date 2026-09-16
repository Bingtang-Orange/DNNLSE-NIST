function [Cplus, Cminus] = assemble_global(defs, x_all, n_all)
% ASSEMBLE_GLOBAL Assembles the global block Cauchy operators C^+ and C^-.
%
% Syntax:
%   [Cplus, Cminus] = assemble_global(defs, x_all, n_all)
%
% Inputs:
%   defs  - Cell array of contour definitions {def1, def2, ..., def_M}
%   x_all - Concatenated collocation nodes on [-1, 1] across all contours
%   n_all - Vector of collocation node counts [n1, n2, ..., nM]
%
% Outputs:
%   Cplus  - Global discrete boundary Cauchy operator matrix C^+ (size N x N)
%   Cminus - Global discrete boundary Cauchy operator matrix C^- (size N x N)
%            satisfying the Plemelj jump relation C^+ - C^- = I on diagonal blocks
%
% Description:
%   Constructs the full multi-contour singular integral operator matrices by
%   iterating through all source-target pairs via assemble_cauchy_block.

    num_curves = length(defs);
    N = sum(n_all);
    Cplus = zeros(N, N);

    block_start = cumsum([0, n_all(1:end-1)]) + 1;
    block_end   = cumsum(n_all);

    % Assemble pairwise interaction blocks
    for i = 1:num_curves
        idx_tgt = block_start(i):block_end(i);
        x_tgt   = x_all(idx_tgt);
        for j = 1:num_curves
            idx_src = block_start(j):block_end(j);
            is_self = (i == j);
            Cplus(idx_tgt, idx_src) = assemble_cauchy_block(x_tgt, n_all(j), defs{j}, defs{i}, is_self);
        end
    end

    % Plemelj jump relation: C^+ - C^- = I (active only on self-interaction diagonal blocks)
    Cminus = Cplus;
    for k = 1:num_curves
        idx = block_start(k):block_end(k);
        Cminus(idx, idx) = Cplus(idx, idx) - eye(n_all(k));
    end
end