function cut = oracle_polar_cut(instance, tri, cell_id, varargin)
%ORACLE_POLAR_CUT Construct and validate one lifted cell-specific polar cut.

parser=inputParser;
addParameter(parser,'Coordinates','source',@(x)ischar(x)||isstring(x));
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(parser,varargin{:});
coordinates=lower(char(parser.Results.Coordinates));
idx=tri.cells(cell_id,:);
d=instance.metadata.n+1;
R0=zeros(d,d); lambda=zeros(d,1); ray_labels=cell(d,1);
for j=1:d-1
    S=tri.labels(idx(j),:);
    [R0(:,j),~]=canonical_lifted_ray(instance,S);
    lambda(j)=lambda_oracle(instance,S);
    ray_labels{j}=S;
end
[R0(:,d),~]=canonical_lifted_ray(instance,'Z');
lambda(d)=lambda_oracle(instance,'Z'); ray_labels{d}='Z';
if strcmp(coordinates,'transformed')
    R=instance.generation.Mx*R0; xbar=instance.pgm_x;
elseif strcmp(coordinates,'source')
    R=R0; xbar=zeros(d,1);
else
    error('oracle_polar_cut:InvalidCoordinates','Coordinates must be source or transformed.');
end
s=svd(R); rank_tol=max(size(R))*eps(norm(R,2));
a=R'\(1./lambda);
lambda_from_cut=1./(a'*R).';
ray_error=abs(lambda_from_cut-lambda)./lambda;

xglobal0=instance.generation.xglobal0;
xalpha0=instance.generation.xalpha0;
if strcmp(coordinates,'transformed')
    xglobal=instance.solution_x; xalpha=instance.incumbent_x;
else
    xglobal=xglobal0; xalpha=xalpha0;
end
mu_global=R\(xglobal-xbar);
mu_alpha=R\(xalpha-xbar);
contains_global=all(mu_global>=-parser.Results.Tolerance);
contains_incumbent=all(mu_alpha>=-parser.Results.Tolerance);
global_lhs=a'*(xglobal-xbar);
incumbent_lhs=a'*(xalpha-xbar);
theory_global=4/(3-instance.metadata.alpha_source);
cut=struct('cell_id',cell_id,'ray_indices',idx,'ray_labels',{ray_labels}, ...
    'R',R,'lambda_oracle',lambda,'a',a,'beta',1, ...
    'lambda_from_cut',lambda_from_cut,'relative_depth_error',ray_error, ...
    'max_relative_depth_error',max(ray_error), ...
    'D_min',min(lambda),'D_mean',mean(lambda), ...
    'D_geo',exp(mean(log(lambda))),'log_volume_multiplier',sum(log(lambda)), ...
    'rank',sum(s>rank_tol),'rank_tolerance',rank_tol,'sigma_min',s(end), ...
    'contains_global',contains_global,'contains_incumbent',contains_incumbent, ...
    'global_lhs',global_lhs,'global_lhs_theory',theory_global, ...
    'global_margin',global_lhs-1,'incumbent_lhs',incumbent_lhs);
cut.passed=cut.rank==d&&cut.max_relative_depth_error<=100*parser.Results.Tolerance&& ...
    (~contains_global||global_lhs>=1-parser.Results.Tolerance)&& ...
    (~contains_incumbent||abs(incumbent_lhs-1)<=100*parser.Results.Tolerance);
end
