function [x, M_prime] = forward_map(def, z)
% FORWARD_MAP Evaluates forward conformal map M(z) and its complex derivative M'(z).
%
% Syntax:
%   [x, M_prime] = forward_map(def, z)
%
% Inputs:
%   def     - Contour definition struct (generated via create_curve_def)
%   z       - Complex coordinate(s) on or near the physical contour (column vector)
%
% Outputs:
%   x       - Conformal coordinate(s) on the reference interval [-1, 1]
%   M_prime - Complex derivative dM/dz evaluated at z
%
% Description:
%   Computes the direct conformal transformation x = M(z) from the physical complex
%   plane back to the unit interval [-1, 1]:
%     * Interval:     M(z) = (a + b - 2*z) / (a - b)
%     * Forward ray:  M(z) = (z - a - c) / (z - a + c),      where c = exp(i*theta)*L
%     * Backward ray: M(z) = (a + c - z) / (z - a + c)
%     * Circular arc: M(z) = K * (z - a - c) / (z - a + c),  where K = -i*cot((th2-th1)/4)

    z = z(:);
    type = lower(def.type);

    switch type
        case 'interval'
            a = def.a;
            b = def.b;
            x = (a + b - 2*z) / (a - b);
            M_prime = (-2 / (a - b)) * ones(size(z));

        case 'ray'
            a     = def.a;
            theta = def.theta;
            L     = def.L;
            c     = exp(1i * theta) * L;

            if ~isfield(def, 'direction') || isempty(def.direction)
                direction = 'forward';
            else
                direction = def.direction;
            end

            denom = z - a + c;
            if strcmp(direction, 'forward')
                x = (z - a - c) ./ denom;
                M_prime = (2 * c) ./ (denom.^2);
            else
                x = (a + c - z) ./ denom;
                M_prime = (-2 * c) ./ (denom.^2);
            end

        case 'arc'
            a   = def.a; 
            r   = def.r;
            th1 = def.theta1; 
            th2 = def.theta2;

            mid_theta = (th1 + th2) / 2;
            c = r * exp(1i * mid_theta);
            coeff = -1i * cot((th2 - th1) / 4);

            denom = z - a + c;
            x = coeff .* (z - a - c) ./ denom;
            M_prime = (2 * coeff * c) ./ (denom.^2);

        otherwise
            error('forward_map: Unsupported contour type "%s".', type);
    end
end