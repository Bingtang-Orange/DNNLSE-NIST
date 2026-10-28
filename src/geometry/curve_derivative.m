function [dz_dx, log_dM_L, log_dM_R] = curve_derivative(def, x)
%CURVE_DERIVATIVE Evaluate conformal mapping derivatives and endpoint metric factors.
%   [DZ_DX, LOG_DM_L, LOG_DM_R] = CURVE_DERIVATIVE(DEF, X) evaluates the
%   derivative dz/dx and endpoint logarithmic metric factors for contour DEF.
%
%   Inputs:
%       def      - Contour definition structure
%       x        - Collocation nodes on the reference interval [-1, 1] (column vector)
%
%   Outputs:
%       dz_dx    - Derivative dz/dx evaluated at x
%       log_dM_L - Logarithmic metric derivative at the left endpoint (x = -1)
%       log_dM_R - Logarithmic metric derivative at the right endpoint (x = +1)

    x = x(:);
    type = lower(def.type);

    switch type
        case 'interval'
            a = def.a;
            b = def.b;
            dz_dx = ((b - a) / 2) * ones(size(x));
            log_dM_L = -log(abs(dz_dx(1)));
            log_dM_R = -log(abs(dz_dx(end)));

        case 'ray'
            a     = def.a;
            theta = def.theta;
            L_ray = def.L;
            c     = exp(1i * theta) * L_ray;

            if isfield(def, 'direction') && ~isempty(def.direction)
                direction = def.direction;
            else
                direction = 'forward';
            end

            if strcmp(direction, 'forward')
                dz_dx = (2 * c) ./ ((1 - x).^2);
                log_dM_L = -log(abs(dz_dx(1)));
                log_dM_R = -Inf;
            else
                dz_dx = (-2 * c) ./ ((1 + x).^2);
                log_dM_L = -Inf;
                log_dM_R = -log(abs(dz_dx(end)));
            end

        case 'arc'
            a   = def.a;
            r   = def.r;
            th1 = def.theta1;
            th2 = def.theta2;

            c = r * exp(1i * (th1 + th2) / 2);
            K = -1i * cot((th2 - th1) / 4);

            dz_dx = (2 * K * c) ./ ((x - K).^2);

            log_dM_L = -log(abs(dz_dx(1)));
            log_dM_R = -log(abs(dz_dx(end)));

        otherwise
            error('curve_derivative: Unsupported contour type "%s".', type);
    end
end