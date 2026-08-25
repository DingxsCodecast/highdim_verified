function instance = generate_highdim_dblp(n, varargin)
%GENERATE_HIGHDIM_DBLP Generate the structural Delta(n,q) DBLP benchmark.
%
% INSTANCE = GENERATE_HIGHDIM_DBLP(N,...) returns the apex-optimal
% structural mode described in the paper. The x-space has dimension N.

validateattributes(n, {'numeric'}, {'scalar','integer','>=',4}, mfilename, 'n');
parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser,'Q',2,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&rem(x,1)==0);
addParameter(parser,'YDimension',n,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
addParameter(parser,'Seed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(parser,'BilinearStrength',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(parser,'Orientation','HH',@(x)ischar(x)||isstring(x));
addParameter(parser,'HouseholderReflections',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
addParameter(parser,'ScalingRange',[0.8 1.2],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>0&&x(2)>=x(1));
addParameter(parser,'NonnegativeMargin',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(parser,'MaterializeXVertices','auto',@(x)islogical(x)||ischar(x)||isstring(x));
addParameter(parser,'MaxVertexElements',5e6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(parser,'Validate',true,@(x)islogical(x)&&isscalar(x));
addParameter(parser,'SaveFile','',@(x)ischar(x)||isstring(x));
parse(parser,varargin{:});
o = parser.Results;
q = o.Q;
if q < 2 || q > n-2
    error('generate_highdim_dblp:InvalidQ','Q must satisfy 2 <= Q <= n-2.');
end
orientation = upper(char(o.Orientation));
ny = o.YDimension;
number_base = nchoosek(n,q);
materialize = resolve_materialization(o.MaterializeXVertices, ...
    n*(number_base+1),o.MaxVertexElements);

old_rng = rng;
cleanup_rng = onCleanup(@()rng(old_rng));
rng(o.Seed,'twister');

[A0,rhs_A0] = build_pyramid_hrep(n,q);
B0 = [-eye(ny);ones(1,ny)];
rhs_B0 = [zeros(ny,1);1];
Vy0 = [zeros(ny,1),eye(ny)];

Qraw = randn(n,ny);
augmented = [Qraw(1:n-1,:);zeros(1,ny)];
[low,high] = qsubset_extrema(augmented,q);
low = low + Qraw(n,:);
high = high + Qraw(n,:);
gamma = max(abs([low,high]),[],'all');
if gamma <= eps
    error('generate_highdim_dblp:ZeroRandomInteraction','Random interaction bound is zero.');
end
Q0 = (o.BilinearStrength/gamma)*Qraw;
scaled_augmented = [Q0(1:n-1,:);zeros(1,ny)];
[interaction_low,interaction_high] = qsubset_extrema(scaled_augmented,q);
interaction_low = interaction_low + Q0(n,:);
interaction_high = interaction_high + Q0(n,:);
interaction_min = min(interaction_low);
interaction_max = max(interaction_high);
interaction_bound = max(abs([interaction_min,interaction_max]));

c0 = [zeros(n-1,1);1];
d0 = ones(ny,1);
certified_gap = min(1,2+interaction_min);
if certified_gap <= 0
    error('generate_highdim_dblp:PlantingFailure','Strict optimality was not preserved.');
end

Mx = make_affine_map(n,orientation,o.HouseholderReflections,o.ScalingRange);
My = make_affine_map(ny,orientation,o.HouseholderReflections,o.ScalingRange);
if strcmp(orientation,'NONE')
    shift_x=zeros(n,1); shift_y=zeros(ny,1);
else
    shift_x = o.NonnegativeMargin-pyramid_row_minima(Mx,n,q);
    shift_y = o.NonnegativeMargin-min([zeros(ny,1),My],[],2);
end

A = A0/Mx;
rhs_A = rhs_A0+A*shift_x;
B = B0/My;
rhs_B = rhs_B0+B*shift_y;
Q = Mx'\(Q0/My);
c = Mx'\c0-Q*shift_y;
d = My'\d0-Q'*shift_x;
solution_x = shift_x;
solution_y = shift_y;
minimum = objective_value(solution_x,solution_y,Q,c,d);

if materialize
    [Vx0,labels] = materialize_pyramid_vertices(n,q);
    Vx = Mx*Vx0+shift_x;
else
    Vx0 = [];
    Vx = [];
    labels = [];
end
Vy = My*Vy0+shift_y;

instance = struct();
instance.mode = 'structural';
instance.Q = Q; instance.A = A; instance.B = B;
instance.rhs_A = rhs_A; instance.rhs_B = rhs_B;
instance.c = c; instance.d = d;
instance.minimum = minimum;
instance.solution_x = solution_x; instance.solution_y = solution_y;
instance.Vx = Vx; instance.Vy = Vy;
instance.metadata = struct( ...
    'generator','generate_highdim_dblp', ...
    'geometry','pyramid_over_Delta(n,q)', ...
    'n',n,'q',q,'ny',ny,'seed',o.Seed, ...
    'orientation',orientation, ...
    'number_x_vertices',number_base+1, ...
    'number_y_vertices',ny+1, ...
    'number_apex_neighbors',number_base, ...
    'number_apex_active_constraints',2*n, ...
    'edge_degeneracy_excess',number_base-n, ...
    'active_constraint_excess',n, ...
    'cross_section_vertex_degree',q*(n-q), ...
    'simple_section_degree',n-1, ...
    'simple_section_deficit',(q-1)*(n-q-1), ...
    'materialized_x_vertices',materialize, ...
    'householder_reflections',o.HouseholderReflections, ...
    'scaling_range',o.ScalingRange, ...
    'nonnegative_margin',o.NonnegativeMargin);
instance.certificate = struct( ...
    'type','analytic q-subset vertex-pair certificate', ...
    'interaction_min',interaction_min, ...
    'interaction_max',interaction_max, ...
    'interaction_bound',interaction_bound, ...
    'certified_vertex_gap',certified_gap, ...
    'unique_global_minimizer',true);
instance.generation = struct( ...
    'Mx',Mx,'My',My,'shift_x',shift_x,'shift_y',shift_y, ...
    'Q0',Q0,'c0',c0,'d0',d0, ...
    'A0',A0,'rhs_A0',rhs_A0,'B0',B0,'rhs_B0',rhs_B0, ...
    'Vx0',Vx0,'Vy0',Vy0,'subset_labels',labels, ...
    'raw_interaction_bound',gamma);

if o.Validate
    instance.validation = validate_highdim_dblp(instance);
    if ~instance.validation.passed
        error('generate_highdim_dblp:ValidationFailed','%s',instance.validation.message);
    end
else
    instance.validation = struct('passed',NaN,'message','Validation was not requested.');
end

save_file = char(o.SaveFile);
if ~isempty(save_file)
    save(save_file,'instance','-v7.3');
end
end

function value = objective_value(x,y,Q,c,d)
value = x'*Q*y+c'*x+d'*y;
end

function materialize = resolve_materialization(request,elements,limit)
if islogical(request)
    materialize = request;
elseif strcmpi(char(request),'auto')
    materialize = elements <= limit;
else
    error('generate_highdim_dblp:InvalidMaterialization', ...
        'MaterializeXVertices must be true, false, or ''auto''.');
end
end
