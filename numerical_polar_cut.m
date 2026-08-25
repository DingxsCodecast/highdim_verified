function cut = numerical_polar_cut(instance, tri, cell_id, depths, varargin)
%NUMERICAL_POLAR_CUT Construct a cell cut from independently computed depths.
%   DEPTHS must come from GENERIC_RAY_DEPTHS; no analytic depth is used here.

parser = inputParser;
addParameter(parser,'Coordinates','source',@(x)ischar(x)||isstring(x));
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(parser,varargin{:});
coordinates = lower(char(parser.Results.Coordinates));
idx = tri.cells(cell_id,:);
d = instance.metadata.n+1;
R0 = zeros(d,d);
for j = 1:d-1
    [R0(:,j),~] = canonical_lifted_ray(instance,tri.labels(idx(j),:));
end
[R0(:,d),~] = canonical_lifted_ray(instance,'Z');
lambda = [depths.lambda(idx);depths.lambda(end)];
if strcmp(coordinates,'transformed')
    R = instance.generation.Mx*R0;
    xbar = instance.pgm_x;
    xglobal = instance.solution_x;
    xalpha = instance.incumbent_x;
elseif strcmp(coordinates,'source')
    R = R0;
    xbar = instance.generation.xbar0;
    xglobal = instance.generation.xglobal0;
    xalpha = instance.generation.xalpha0;
else
    error('numerical_polar_cut:InvalidCoordinates', ...
        'Coordinates must be source or transformed.');
end
s = svd(R);
rank_tol = max(size(R))*eps(max(1,norm(R,2)));
a = R'\(1./lambda);
lambda_from_cut = 1./(a'*R).';
cut_recovery_error = abs(lambda_from_cut-lambda)./lambda;
mu_global = R\(xglobal-xbar);
mu_alpha = R\(xalpha-xbar);
contains_global = all(mu_global>=-parser.Results.Tolerance);
contains_incumbent = all(mu_alpha>=-parser.Results.Tolerance);
global_lhs = a'*(xglobal-xbar);
incumbent_lhs = a'*(xalpha-xbar);
cut = struct('cell_id',cell_id,'R',R,'lambda_numerical',lambda,'a',a, ...
    'lambda_from_cut',lambda_from_cut,'cut_recovery_error',cut_recovery_error, ...
    'max_cut_recovery_error',max(cut_recovery_error), ...
    'D_min',min(lambda),'D_mean',mean(lambda), ...
    'D_geo',exp(mean(log(lambda))),'log_volume_multiplier',sum(log(lambda)), ...
    'rank',sum(s>rank_tol),'rank_tolerance',rank_tol,'sigma_min',s(end), ...
    'contains_global',contains_global,'contains_incumbent',contains_incumbent, ...
    'global_lhs',global_lhs,'global_margin',global_lhs-1, ...
    'incumbent_lhs',incumbent_lhs);
cut.passed = cut.rank==d&&cut.max_cut_recovery_error<=100*parser.Results.Tolerance&& ...
    (~contains_global||global_lhs>=1-parser.Results.Tolerance)&& ...
    (~contains_incumbent||abs(incumbent_lhs-1)<=100*parser.Results.Tolerance);
end
