function [dz_dx, log_dM_L, log_dM_R] = curve_derivative(def, x)
% CURVE_DERIVATIVE Evaluates conformal mapping derivatives and endpoint metric factors.
%
% Syntax:
%   [dz_dx, log_dM_L, log_dM_R] = curve_derivative(def, x)
%
% Inputs:
%   def      - Contour definition struct (generated via create_curve_def)
%   x        - Collocation nodes on the reference interval [-1, 1] (column vector)
%
% Outputs:
%   dz_dx    - Derivative dz/dx of the inverse map z(x), evaluated at x
%   log_dM_L - log|M'(z_L)| at the left endpoint (x = -1)
%   log_dM_R - log|M'(z_R)| at the right endpoint (x = +1)
%
% Description:
%   Computes metric derivatives required for diagonal endpoint regularization
%   in discrete Cauchy matrices, strictly adhering to Olver (2012), Eq. (4.3).
%   If an endpoint lies at infinity (unbounded rays), the metric derivative is -Inf.

    x = x(:);
    type = lower(def.type);

    switch type
        case 'interval'
            a = def.a;
            b = def.b;
            % Inverse map: z(x) = a + (b-a)*(x+1)/2
            dz_dx = ((b - a) / 2) * ones(size(x));
            % Since M'(z) = 1 / (dz/dx), log|M'| = -log|dz/dx|
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
                % Forward ray: z(x) = a + c*(1+x)/(1-x)
                % dz/dx = 2c / (1-x)^2
                dz_dx = (2 * c) ./ ((1 - x).^2);
                % Left endpoint (x = -1) is finite; right endpoint (x = 1) is at infinity
                log_dM_L = -log(abs(dz_dx(1))); % Evaluates to log(2 / |c|)
                log_dM_R = -Inf;
            else
                % Backward ray: z(x) = a + c*(1-x)/(1+x)
                % dz/dx = -2c / (1+x)^2
                dz_dx = (-2 * c) ./ ((1 + x).^2);
                % Left endpoint (x = -1) is at infinity; right endpoint (x = 1) is finite
                log_dM_L = -Inf;
                log_dM_R = -log(abs(dz_dx(end))); % Evaluates to log(2 / |c|)
            end

        case 'arc'
            a   = def.a;
            r   = def.r;
            th1 = def.theta1;
            th2 = def.theta2;

            % Conformal parameters for circular arc Möbius map
            c = r * exp(1i * (th1 + th2) / 2);
            K = -1i * cot((th2 - th1) / 4);

            % Inverse mapping derivative: dz/dx = 2 * K * c / (x - K)^2
            dz_dx = (2 * K * c) ./ ((x - K).^2);

            % Endpoint logarithmic metric derivatives
            log_dM_L = -log(abs(dz_dx(1)));
            log_dM_R = -log(abs(dz_dx(end)));

        otherwise
            error('curve_derivative: Unsupported contour type "%s".', type);
    end
end