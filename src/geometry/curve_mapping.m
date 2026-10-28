function [curves_points, curves_info] = curve_mapping(curve_defs, x_interval)
%CURVE_MAPPING Map reference interval collocation nodes to physical complex contours.
%   [CURVES_POINTS, CURVES_INFO] = CURVE_MAPPING(CURVE_DEFS, X_INTERVAL)
%   evaluates the conformal inverse map z = M^{-1}(x) from [-1, 1] onto physical
%   segments, circular arcs, or unbounded rays.
%
%   Inputs:
%       curve_defs    - Array or struct of contour definitions
%       x_interval    - Reference collocation nodes on [-1, 1] (column vector)
%
%   Outputs:
%       curves_points - Cell array of mapped complex coordinates for each contour
%       curves_info   - Struct array containing contour metadata and image of infinity

    if nargin < 2
        error('curve_mapping: Requires contour definitions and reference nodes.');
    end

    x = x_interval(:);
    num_curves = length(curve_defs);
    curves_points = cell(num_curves, 1);
    curves_info = struct('M_inf', cell(num_curves, 1), 'type', cell(num_curves, 1));

    for idx = 1:num_curves
        def  = curve_defs(idx);
        type = lower(def.type);
        curves_info(idx).type = type;

        switch type
            case 'interval'
                a = def.a;
                b = def.b;
                z = a + (b - a) * (x + 1) / 2;
                curves_points{idx} = z;

                if isfield(def, 'M_inf') && ~isempty(def.M_inf)
                    curves_info(idx).M_inf = def.M_inf;
                else
                    curves_info(idx).M_inf = inf;
                end

            case 'ray'
                a     = def.a;
                theta = def.theta;
                if ~isfield(def, 'L') || isempty(def.L)
                    error('curve_mapping: Ray definition must contain scaling parameter L.');
                end
                L_ray = def.L;
                c     = exp(1i * theta) * L_ray;

                if isfield(def, 'direction') && ~isempty(def.direction)
                    direction = def.direction;
                else
                    direction = 'forward';
                end

                if strcmp(direction, 'forward')
                    x_use = x(1:end-1);
                    if isempty(x_use)
                        warning('curve_mapping: All nodes removed for forward ray %d.', idx);
                        curves_points{idx} = [];
                        curves_info(idx).M_inf = 1;
                        continue;
                    end
                    z = a + c * (1 + x_use) ./ (1 - x_use);
                    curves_info(idx).M_inf = 1;
                else
                    x_use = x(2:end);
                    if isempty(x_use)
                        warning('curve_mapping: All nodes removed for backward ray %d.', idx);
                        curves_points{idx} = [];
                        curves_info(idx).M_inf = -1;
                        continue;
                    end
                    z = a + c * (1 - x_use) ./ (1 + x_use);
                    curves_info(idx).M_inf = -1;
                end
                curves_points{idx} = z;

            case 'arc'
                a   = def.a; 
                r   = def.r;
                th1 = def.theta1; 
                th2 = def.theta2;
                x_use = x(:);

                mid_theta = (th1 + th2) / 2;
                c = r * exp(1i * mid_theta);
                coeff = -1i * cot((th2 - th1) / 4);

                z = a + c .* (x_use + coeff) ./ (coeff - x_use);

                curves_points{idx} = z;
                curves_info(idx).M_inf = coeff;

            otherwise
                error('curve_mapping: Unsupported contour type "%s".', type);
        end
    end
end