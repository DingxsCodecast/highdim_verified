function depths = generic_ray_depths(instance, labels, varargin)
%GENERIC_RAY_DEPTHS Compute all subset-ray and z-ray depths by black-box bisection.

timer = tic;
nr = size(labels,1);
lambda = zeros(nr+1,1);
residual = zeros(nr+1,1);
iterations = zeros(nr+1,1);
evaluations = zeros(nr+1,1);
bracket_steps = zeros(nr+1,1);
for k = 1:nr
    root = generic_ray_extension(instance,labels(k,:),varargin{:});
    lambda(k) = root.lambda;
    residual(k) = root.residual;
    iterations(k) = root.iterations;
    evaluations(k) = root.function_evaluations;
    bracket_steps(k) = root.bracket_steps;
end
root = generic_ray_extension(instance,'Z',varargin{:});
lambda(end) = root.lambda;
residual(end) = root.residual;
iterations(end) = root.iterations;
evaluations(end) = root.function_evaluations;
bracket_steps(end) = root.bracket_steps;
depths = struct('labels',labels,'lambda',lambda,'residual',residual, ...
    'iterations',iterations,'function_evaluations',evaluations, ...
    'bracket_steps',bracket_steps,'seconds',toc(timer), ...
    'passed',all(isfinite(lambda))&&all(lambda>1));
end
