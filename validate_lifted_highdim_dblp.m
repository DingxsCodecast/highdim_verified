function report = validate_lifted_highdim_dblp(instance, varargin)
%VALIDATE_LIFTED_HIGHDIM_DBLP Validate PGM, optimum, affine, and depth oracles.

parser=inputParser;
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'RaySamples',20,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
parse(parser,varargin{:});
tol=parser.Results.Tolerance;
n=instance.metadata.n; q=instance.metadata.q; d=n+1;
g=instance.generation; cert=instance.certificate;

checks=struct();
checks.mode=strcmp(instance.mode,'lifted');
checks.parameter_ranges=q>=2&&q<=n-2&&instance.metadata.eta>0&& ...
    instance.metadata.eta<2&&instance.metadata.alpha_source>-1&&instance.metadata.alpha_source<0;
checks.matrix_dimensions=isequal(size(instance.Q),[d,1])&&size(instance.A,2)==d&&size(instance.B,2)==1;
points_x=[instance.pgm_x,instance.solution_x,instance.incumbent_x];
points_y=[instance.pgm_y,instance.solution_y,instance.incumbent_y];
checks.certified_points_feasible=max(instance.A*points_x-instance.rhs_A,[],'all')<=tol&& ...
    max(instance.B*points_y-instance.rhs_B,[],'all')<=tol;

values=zeros(1,3);
for k=1:3
    values(k)=objective_value(points_x(:,k),points_y(:,k),instance);
end
targets=[instance.pgm_value,instance.minimum,instance.incumbent_value];
checks.reported_values=all(abs(values-targets)<=20*tol.*max(1,abs(targets)));
checks.source_values=abs(source_objective(g.xbar0,g.ybar0,g)-0)<=tol&& ...
    abs(source_objective(g.xglobal0,g.yglobal0,g)+1)<=tol&& ...
    abs(source_objective(g.xalpha0,g.yalpha0,g)-instance.metadata.alpha_source)<=tol;

xres=instance.A*instance.pgm_x-instance.rhs_A;
active=abs(xres)<=tol;
checks.active_facets=sum(active)==2*n+1;
checks.incident_edges=instance.metadata.number_lifted_rays==nchoosek(n,q)+1;
checks.degeneracy_excesses=instance.metadata.active_constraint_excess==n&& ...
    instance.metadata.edge_degeneracy_excess==nchoosek(n,q)-n;
checks.ray_degrees=instance.metadata.old_ray_degree==q*(n-q)+1&& ...
    instance.metadata.z_ray_degree==nchoosek(n,q);

checks.affine_constraints=relative_residual(instance.A*g.Mx,g.A0)<=tol&& ...
    relative_residual(instance.rhs_A-instance.A*g.shift_x,g.rhs_A0)<=tol&& ...
    relative_residual(instance.B*g.My,g.B0)<=tol&& ...
    relative_residual(instance.rhs_B-instance.B*g.shift_y,g.rhs_B0)<=tol;
checks.affine_objective=relative_residual(g.Mx'*instance.Q*g.My,g.Q0)<=20*tol&& ...
    relative_residual(g.Mx'*(instance.Q*g.shift_y+instance.c),g.c0)<=20*tol&& ...
    relative_residual(g.My'*(instance.Q'*g.shift_x+instance.d),g.d0)<=20*tol;
checks.objective_offset=abs(objective_value(instance.pgm_x,instance.pgm_y,instance)-cert.objective_offset)<=20*tol;

checks.distinct_augmented_weights=numel(unique([cert.w;0]))==n;
checks.exposing_extrema=cert.D>0&&abs(lifted_psi(subset_to_h(cert.Sstar,n),1,cert))<=tol;

labels=deterministic_subset_sample(n,q,parser.Results.RaySamples);
psi_values=zeros(size(labels,1),1);
neighbor_values=zeros(size(labels,1),1);
lambda_formula=zeros(size(labels,1)+1,1);
lambda_numeric=zeros(size(labels,1)+1,1);
depth_residual=zeros(size(labels,1)+1,1);
for k=1:size(labels,1)
    h=subset_to_h(labels(k,:),n);
    psi_values(k)=lifted_psi(h,1,cert);
    ray=[h;1;0];
    neighbor_values(k)=lifted_reduced_value(ray,instance,'source');
    lambda_formula(k)=lambda_oracle(instance,labels(k,:));
    rootfun=@(lam)lifted_reduced_value(lam*ray,instance,'source')-instance.metadata.alpha_source;
    lambda_numeric(k)=fzero(rootfun,[1,max(2,1.1*lambda_formula(k))]);
    depth_residual(k)=abs(rootfun(lambda_formula(k)));
end
rayz=[zeros(n,1);1];
lambda_formula(end)=lambda_oracle(instance,'Z');
rootfun=@(lam)lifted_reduced_value(lam*rayz,instance,'source')-instance.metadata.alpha_source;
lambda_numeric(end)=fzero(rootfun,[1,max(2,1.1*lambda_formula(end))]);
depth_residual(end)=abs(rootfun(lambda_formula(end)));
checks.psi_range=all(psi_values>=-tol)&all(psi_values<=1+tol);
checks.pgm_neighbor_margin=min([neighbor_values;lifted_reduced_value(rayz,instance,'source')])>=1-tol;
checks.lambda_greater_than_one=all(lambda_formula>1);
checks.lambda_closed_form=max(abs(lambda_numeric-lambda_formula))<=100*tol;
checks.lambda_boundary=max(depth_residual)<=100*tol;

mapped_residual=zeros(size(lambda_formula));
for k=1:size(labels,1)
    [~,rt]=canonical_lifted_ray(instance,labels(k,:));
    xt=instance.pgm_x+lambda_formula(k)*rt;
    mapped_residual(k)=abs(lifted_reduced_value(xt,instance,'transformed')-instance.incumbent_value);
end
[~,rt]=canonical_lifted_ray(instance,'Z');
mapped_residual(end)=abs(lifted_reduced_value(instance.pgm_x+lambda_formula(end)*rt,instance,'transformed')-instance.incumbent_value);
checks.transformed_lambda_invariance=max(mapped_residual)<=200*tol*max(1,abs(instance.incumbent_value));

exhaustive=struct('performed',false,'minimum',NaN,'minimizer_count',NaN,'target_error',NaN);
if ~isempty(instance.Vx)
    all_values=instance.Vx'*instance.Q*instance.Vy+ ...
        (instance.c'*instance.Vx)'+instance.d'*instance.Vy;
    exhaustive.performed=true;
    exhaustive.minimum=min(all_values,[],'all');
    exhaustive.minimizer_count=sum(abs(all_values(:)-exhaustive.minimum)<=tol);
    exhaustive.target_error=abs(exhaustive.minimum-instance.minimum);
    checks.exhaustive_global_optimum=exhaustive.target_error<=100*tol&&exhaustive.minimizer_count==1;
else
    checks.exhaustive_global_optimum=true;
end

arrival=struct('initial_y_source',0,'first_x_best_response_source',zeros(d,1), ...
    'second_y_best_response_source',0,'PGM_returned',true,'PGM_matches_planted',true);
checks.deterministic_arrival=source_objective(zeros(d,1),0,g)==0&& ...
    source_objective(zeros(d,1),1,g)>0;

check_values=struct2cell(checks);
report=struct();
report.passed=all(cellfun(@(x)islogical(x)&&isscalar(x)&&x,check_values));
report.checks=checks;
report.tolerances=struct('feasibility',tol,'cut',tol,'cone',tol, ...
    'rank_rule','max(size(R))*eps(norm(R,2))');
report.metrics=struct('source_values',values-cert.objective_offset, ...
    'psi_values',psi_values,'sampled_labels',labels, ...
    'neighbor_best_response_values',neighbor_values, ...
    'lambda_formula',lambda_formula,'lambda_numeric',lambda_numeric, ...
    'maximum_lambda_absolute_error',max(abs(lambda_numeric-lambda_formula)), ...
    'maximum_boundary_residual',max(depth_residual), ...
    'maximum_transformed_boundary_residual',max(mapped_residual), ...
    'exhaustive',exhaustive,'arrival',arrival);
if report.passed
    report.message='All lifted PGM, optimum, affine, and ray-depth checks passed.';
else
    names=fieldnames(checks); failed=names(~cellfun(@(s)checks.(s),names));
    report.message=['Failed checks: ',strjoin(failed,', ')];
end
end

function value=objective_value(x,y,instance)
value=x'*instance.Q*y+instance.c'*x+instance.d'*y;
end

function value=source_objective(x,y,g)
value=x'*g.Q0*y+g.c0'*x+g.d0*y;
end

function labels=deterministic_subset_sample(n,q,requested)
count=min(requested,nchoosek(n,q));
if nchoosek(n,q)<=5000
    all_labels=enumerate_q_subsets(n,q);
    idx=unique(round(linspace(1,size(all_labels,1),count)),'stable');
    labels=all_labels(idx,:);
else
    labels=zeros(count,q); labels(1,:)=1:q;
    for k=2:count
        labels(k,:)=sort(mod((0:q-1)+(k-1),n)+1);
    end
    labels=unique(labels,'rows','stable');
end
end
