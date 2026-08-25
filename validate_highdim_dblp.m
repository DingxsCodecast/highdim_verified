function report = validate_highdim_dblp(instance, varargin)
%VALIDATE_HIGHDIM_DBLP Validate the structural Delta(n,q) benchmark.

parser = inputParser;
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(parser,'EdgeSamples',20,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
parse(parser,varargin{:});
tol = parser.Results.Tolerance;
n = instance.metadata.n;
q = instance.metadata.q;
ny = instance.metadata.ny;

checks = struct();
checks.mode = strcmp(instance.mode,'structural');
checks.matrix_dimensions = isequal(size(instance.Q),[n,ny]) && ...
    size(instance.A,2)==n && size(instance.B,2)==ny;
x_residual = instance.A*instance.solution_x-instance.rhs_A;
y_residual = instance.B*instance.solution_y-instance.rhs_B;
checks.planted_solution_feasible = max(x_residual)<=tol && max(y_residual)<=tol;
computed = objective_value(instance.solution_x,instance.solution_y,instance);
checks.reported_value_consistent = abs(computed-instance.minimum)<=tol*max(1,abs(instance.minimum));

active = abs(x_residual)<=tol;
active_count = sum(active);
checks.active_facet_certificate = active_count==2*n;
checks.incident_edge_certificate = instance.metadata.number_apex_neighbors==nchoosek(n,q);
checks.degeneracy_excesses = instance.metadata.active_constraint_excess==n && ...
    instance.metadata.edge_degeneracy_excess==nchoosek(n,q)-n;
checks.section_degree = instance.metadata.cross_section_vertex_degree==q*(n-q) && ...
    instance.metadata.simple_section_deficit==(q-1)*(n-q-1);

g = instance.generation;
checks.affine_constraints = relative_residual(instance.A*g.Mx,g.A0)<=tol && ...
    relative_residual(instance.rhs_A-instance.A*g.shift_x,g.rhs_A0)<=tol && ...
    relative_residual(instance.B*g.My,g.B0)<=tol && ...
    relative_residual(instance.rhs_B-instance.B*g.shift_y,g.rhs_B0)<=tol;
checks.affine_objective = relative_residual(g.Mx'*instance.Q*g.My,g.Q0)<=20*tol && ...
    relative_residual(g.Mx'*(instance.Q*g.shift_y+instance.c),g.c0)<=20*tol && ...
    relative_residual(g.My'*(instance.Q'*g.shift_x+instance.d),g.d0)<=20*tol;

implicit_x_min = pyramid_row_minima(g.Mx,n,q)+g.shift_x;
implicit_y_min = min([zeros(ny,1),g.My],[],2)+g.shift_y;
checks.nonnegative_vertices = min(implicit_x_min)>=-tol && min(implicit_y_min)>=-tol;

basis_labels = independent_ray_labels(n,q);
basis = zeros(n,n);
for k=1:n
    basis(:,k) = [subset_to_h(basis_labels(k,:),n);1];
end
rank_tol = max(size(basis))*eps(norm(basis,2));
checks.full_dimensional = sum(svd(basis)>rank_tol)==n;

sample_labels = deterministic_subset_sample(n,q,parser.Results.EdgeSamples);
edge_ranks = zeros(size(sample_labels,1),1);
neighbor_counts = zeros(size(sample_labels,1),1);
for k=1:size(sample_labels,1)
    x0 = [subset_to_h(sample_labels(k,:),n);1];
    x = g.Mx*x0+g.shift_x;
    common = active & abs(instance.A*x-instance.rhs_A)<=tol;
    s = svd(instance.A(common,:),'econ');
    local_tol = max(size(instance.A(common,:)))*eps(max(1,s(1)));
    edge_ranks(k) = sum(s>local_tol);
    neighbor_counts(k) = size(delta_neighbors(sample_labels(k,:),n),1);
end
checks.sampled_incident_edges = all(edge_ranks==n-1);
checks.sampled_neighborhoods = all(neighbor_counts==q*(n-q));

augmented = [g.Q0(1:n-1,:);zeros(1,ny)];
[low,high] = qsubset_extrema(augmented,q);
low = low+g.Q0(n,:); high = high+g.Q0(n,:);
interaction_min = min(low); interaction_max = max(high);
interaction_bound = max(abs([interaction_min,interaction_max]));
gap = min(1,2+interaction_min);
checks.interaction_certificate = interaction_bound<1 && ...
    abs(interaction_min-instance.certificate.interaction_min)<=tol && ...
    abs(interaction_max-instance.certificate.interaction_max)<=tol;
checks.global_gap_certificate = gap>0 && ...
    abs(gap-instance.certificate.certified_vertex_gap)<=tol && ...
    instance.certificate.unique_global_minimizer;

exhaustive = struct('performed',false,'minimum',NaN,'gap',NaN, ...
    'minimizer_count',NaN,'target_error',NaN);
if ~isempty(instance.Vx)
    values = instance.Vx'*instance.Q*instance.Vy + ...
        (instance.c'*instance.Vx)' + instance.d'*instance.Vy;
    target = values(1,1);
    others = values(:); others(1)=[];
    exhaustive.performed = true;
    exhaustive.minimum = min(values,[],'all');
    exhaustive.gap = min(others)-target;
    exhaustive.minimizer_count = sum(abs(values(:)-exhaustive.minimum)<=tol);
    exhaustive.target_error = abs(target-instance.minimum);
    checks.exhaustive_vertex_pairs = exhaustive.target_error<=50*tol && ...
        exhaustive.gap>=instance.certificate.certified_vertex_gap-50*tol && ...
        exhaustive.minimizer_count==1;
    checks.materialized_feasibility = ...
        max(instance.A*instance.Vx-instance.rhs_A,[],'all')<=tol && ...
        max(instance.B*instance.Vy-instance.rhs_B,[],'all')<=tol;
else
    checks.exhaustive_vertex_pairs = true;
    checks.materialized_feasibility = true;
end

values = struct2cell(checks);
report = struct();
report.passed = all(cellfun(@(x)islogical(x)&&isscalar(x)&&x,values));
report.checks = checks;
report.tolerances = struct('feasibility',tol,'algebra',tol,'rank_rule','max(size(R))*eps(norm(R,2))');
report.metrics = struct('maximum_x_violation',max(x_residual), ...
    'maximum_y_violation',max(y_residual), ...
    'active_constraints_at_apex',active_count, ...
    'sampled_edge_ranks',edge_ranks, ...
    'sampled_neighbor_counts',neighbor_counts, ...
    'sampled_labels',sample_labels, ...
    'exhaustive',exhaustive);
if report.passed
    report.message = 'All structural Delta(n,q) checks passed.';
else
    names = fieldnames(checks);
    failed = names(~cellfun(@(s)checks.(s),names));
    report.message = ['Failed checks: ',strjoin(failed,', ')];
end
end

function value = objective_value(x,y,instance)
value=x'*instance.Q*y+instance.c'*x+instance.d'*y;
end

function labels = independent_ray_labels(n,q)
S0 = 1:q;
labels = zeros(n,q);
labels(1,:) = S0;
cursor = 1;
for j=q+1:n
    cursor=cursor+1; T=S0; T(1)=j; labels(cursor,:)=sort(T);
end
for i=2:q
    cursor=cursor+1; T=S0; T(i)=q+1; labels(cursor,:)=sort(T);
end
end

function labels = deterministic_subset_sample(n,q,requested)
count = min(requested,nchoosek(n,q));
all_count = nchoosek(n,q);
if all_count<=5000
    all_labels = enumerate_q_subsets(n,q);
    index = unique(round(linspace(1,all_count,count)),'stable');
    labels = all_labels(index,:);
else
    labels = zeros(count,q);
    labels(1,:) = 1:q;
    for k=2:count
        shift = mod(k-1,n);
        labels(k,:) = sort(mod((0:q-1)+shift,n)+1);
    end
    labels = unique(labels,'rows','stable');
end
end
