function def = create_curve_def(type, varargin)
% CREATE_CURVE_DEF Initializes a geometric contour specification struct.
%
% Syntax:
%   def = create_curve_def('interval', a, b)
%   def = create_curve_def('ray', a, theta, L)
%   def = create_curve_def('ray', a, theta, L, direction) % 'forward' (default) or 'backward'
%   def = create_curve_def('arc', a, r, theta1, theta2)
%
% Inputs:
%   type      - String specifying contour family: 'interval', 'ray', or 'arc'
%   varargin  - Parameter list corresponding to the chosen contour type:
%                 * 'interval': a (start point), b (end point)
%                 * 'ray'     : a (origin), theta (launch angle), L (scaling parameter),
%                               direction ('forward' or 'backward')
%                 * 'arc'     : a (center), r (radius), theta1 (start angle), theta2 (end angle)
%
% Outputs:
%   def       - Standardized contour definition struct
%
% References:
%   S. Olver, Numer. Math. (2012), Section 4.
%   T. Trogdon and S. Olver, "Riemann-Hilbert Problems, Their Numerical Solution,
%   and the Computation of Nonlinear Special Functions", SIAM (2016).

    % Initialize all struct fields to maintain consistent memory layout
    def = struct('type', type, ...
                 'a', [], 'b', [], ...
                 'theta', [], 'L', [], 'direction', [], ...
                 'r', [], 'theta1', [], 'theta2', [], ...
                 'M_inf', []);

    switch lower(type)
        case 'interval'
            if nargin ~= 3
                error('create_curve_def: Interval requires exactly two endpoints: a and b.');
            end
            def.a = varargin{1};
            def.b = varargin{2};

        case 'ray'
            if nargin < 4
                error('create_curve_def: Ray requires at least origin (a), angle (theta), and scale (L).');
            end
            def.a     = varargin{1};
            def.theta = varargin{2};
            def.L     = varargin{3};
            
            if nargin >= 5 && ~isempty(varargin{4})
                def.direction = varargin{4};
            else
                def.direction = 'forward'; % Default orientation
            end

        case 'arc'
            if nargin ~= 5
                error('create_curve_def: Arc requires center (a), radius (r), theta1, and theta2.');
            end
            def.a      = varargin{1};
            def.r      = varargin{2};
            def.theta1 = varargin{3};
            def.theta2 = varargin{4};

        otherwise
            error('create_curve_def: Unsupported contour type "%s".', type);
    end
end