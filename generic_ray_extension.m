function result = generic_ray_extension(instance, ray_label, varargin)
%GENERIC_RAY_EXTENSION Black-box bisection for g(xbar+lambda*r)=alpha.
%   This routine deliberately does not call LAMBDA_ORACLE or use its formula.

parser = inputParser;
addParameter(parser,'Coordinates','source',@(x)ischar(x)||isstring(x));
addParameter(parser,'RelativeTolerance',1e-13,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'FunctionTolerance',1e-13,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'InitialUpper',1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'MaxBracketSteps',60,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
addParameter(parser,'MaxIterations',200,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
parse(parser,varargin{:});
o = parser.Results;
coordinates = lower(char(o.Coordinates));

[ray0,rayt] = canonical_lifted_ray(instance,ray_label);
if strcmp(coordinates,'source')
    xbar = instance.generation.xbar0;
    ray = ray0;
    alpha = instance.metadata.alpha_source;
elseif strcmp(coordinates,'transformed')
    xbar = instance.pgm_x;
    ray = rayt;
    alpha = instance.incumbent_value;
else
    error('generic_ray_extension:InvalidCoordinates', ...
        'Coordinates must be source or transformed.');
end

evaluate = @(lambda) lifted_reduced_value(xbar+lambda*ray,instance,coordinates)-alpha;
lo = 0;
f_lower = evaluate(lo);
evaluations = 1;
if ~(f_lower > 0)
    error('generic_ray_extension:InvalidInitialSign', ...
        'The planted PGM must lie strictly above the incumbent level.');
end
hi = o.InitialUpper;
f_upper = evaluate(hi);
evaluations = evaluations+1;
bracket_steps = 0;
while f_upper > 0 && bracket_steps < o.MaxBracketSteps
    hi = 2*hi;
    f_upper = evaluate(hi);
    evaluations = evaluations+1;
    bracket_steps = bracket_steps+1;
end
if f_upper > 0
    error('generic_ray_extension:NoBracket', ...
        'Failed to bracket a positive boundary intersection.');
end

iterations = 0;
while iterations < o.MaxIterations
    midpoint = lo+(hi-lo)/2;
    f_midpoint = evaluate(midpoint);
    evaluations = evaluations+1;
    iterations = iterations+1;
    if abs(f_midpoint) <= o.FunctionTolerance
        lo = midpoint;
        hi = midpoint;
        f_lower = f_midpoint;
        f_upper = f_midpoint;
        break;
    elseif f_midpoint > 0
        lo = midpoint;
        f_lower = f_midpoint;
    else
        hi = midpoint;
        f_upper = f_midpoint;
    end
    width = hi-lo;
    if width <= o.RelativeTolerance*max(1,abs(midpoint))
        break;
    end
end
lambda = lo+(hi-lo)/2;
residual = evaluate(lambda);
evaluations = evaluations+1;
result = struct('lambda',lambda,'residual',residual, ...
    'lower',lo,'upper',hi,'f_lower',f_lower,'f_upper',f_upper, ...
    'bracket_steps',bracket_steps,'iterations',iterations, ...
    'function_evaluations',evaluations,'coordinates',coordinates, ...
    'relative_tolerance',o.RelativeTolerance,'function_tolerance',o.FunctionTolerance, ...
    'passed',f_lower>=-o.FunctionTolerance&&f_upper<=o.FunctionTolerance&& ...
        hi-lo<=10*o.RelativeTolerance*max(1,abs(lambda)));
end
