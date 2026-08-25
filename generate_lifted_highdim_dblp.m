function instance = generate_lifted_highdim_dblp(n, varargin)
%GENERATE_LIFTED_HIGHDIM_DBLP Generate the strict non-global-PGM benchmark.
%
% N is the dimension of the underlying pyramid; the returned x-dimension is
% N+1. Canonical rays are adjacent-vertex differences and are never normalized.

validateattributes(n,{'numeric'},{'scalar','integer','>=',4},mfilename,'n');
parser=inputParser;
addParameter(parser,'Q',2,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&rem(x,1)==0);
addParameter(parser,'Seed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(parser,'Eta',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<2);
addParameter(parser,'Alpha',-0.5,@(x)isnumeric(x)&&isscalar(x)&&x>-1&&x<0);
addParameter(parser,'Orientation','HH',@(x)ischar(x)||isstring(x));
addParameter(parser,'HouseholderReflections',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
addParameter(parser,'ScalingRange',[0.8 1.2],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>0&&x(2)>=x(1));
addParameter(parser,'NonnegativeMargin',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(parser,'MaterializeXVertices','auto',@(x)islogical(x)||ischar(x)||isstring(x));
addParameter(parser,'MaxVertexElements',5e6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(parser,'Validate',true,@(x)islogical(x)&&isscalar(x));
addParameter(parser,'SaveFile','',@(x)ischar(x)||isstring(x));
parse(parser,varargin{:});
o=parser.Results;
q=o.Q;
if q<2||q>n-2
    error('generate_lifted_highdim_dblp:InvalidQ','Q must satisfy 2 <= Q <= n-2.');
end
orientation=upper(char(o.Orientation));
number_base=nchoosek(n,q);
materialize=resolve_materialization(o.MaterializeXVertices, ...
    (n+1)*2*(number_base+1),o.MaxVertexElements);

old_rng=rng;
cleanup_rng=onCleanup(@()rng(old_rng));
rng(o.Seed,'twister');

[A0,rhs_A0]=build_lifted_hrep(n,q);
B0=[-1;1]; rhs_B0=[0;1]; Vy0=[0,1];

while true
    w=randn(n-1,1);
    waug=[w;0];
    if numel(unique(waug))==n
        break;
    end
end
[ordered,index]=sort(waug,'ascend');
Sstar=sort(index(1:q)).';
m=sum(ordered(1:q));
M=sum(ordered(end-q+1:end));
D=M-m;
if ~(D>0)
    error('generate_lifted_highdim_dblp:InvalidExposingFunctional','M-m must be positive.');
end

Q0=[zeros(n-1,1);-3;-3];
c0=[o.Eta*w/D;1-o.Eta*m/D;1];
d0=3;
xbar0=zeros(n+1,1); ybar0=0;
hstar=subset_to_h(Sstar,n);
xglobal0=[hstar;1;1]; yglobal0=1;
zalpha=(1-o.Alpha)/2;
xalpha0=[hstar;1;zalpha]; yalpha0=1;

Mx=make_affine_map(n+1,orientation,o.HouseholderReflections,o.ScalingRange);
My=make_affine_map(1,orientation,o.HouseholderReflections,o.ScalingRange);
if strcmp(orientation,'NONE')
    shift_x=zeros(n+1,1); shift_y=0;
else
    shift_x=o.NonnegativeMargin-lifted_row_minima(Mx,n,q);
    shift_y=o.NonnegativeMargin-min([0,My],[],2);
end

A=A0/Mx; rhs_A=rhs_A0+A*shift_x;
B=B0/My; rhs_B=rhs_B0+B*shift_y;
Q=Mx'\(Q0/My);
c=Mx'\c0-Q*shift_y;
d=My'\d0-Q'*shift_x;
kappa=shift_x'*Q*shift_y+c'*shift_x+d'*shift_y;
xbar=Mx*xbar0+shift_x; ybar=My*ybar0+shift_y;
xglobal=Mx*xglobal0+shift_x; yglobal=My*yglobal0+shift_y;
xalpha=Mx*xalpha0+shift_x; yalpha=My*yalpha0+shift_y;

if materialize
    [Vpyr0,labels]=materialize_pyramid_vertices(n,q);
    Vx0=[[Vpyr0;zeros(1,size(Vpyr0,2))], ...
        [Vpyr0;ones(1,size(Vpyr0,2))]];
    Vx=Mx*Vx0+shift_x;
else
    Vx0=[]; Vx=[]; labels=[];
end
Vy=My*Vy0+shift_y;

instance=struct();
instance.mode='lifted';
instance.Q=Q; instance.A=A; instance.B=B;
instance.rhs_A=rhs_A; instance.rhs_B=rhs_B;
instance.c=c; instance.d=d;
instance.pgm_x=xbar; instance.pgm_y=ybar; instance.pgm_value=kappa;
instance.solution_x=xglobal; instance.solution_y=yglobal; instance.minimum=kappa-1;
instance.incumbent_x=xalpha; instance.incumbent_y=yalpha; instance.incumbent_value=kappa+o.Alpha;
instance.Vx=Vx; instance.Vy=Vy;
instance.metadata=struct( ...
    'generator','generate_lifted_highdim_dblp', ...
    'geometry','pyramid_over_Delta(n,q)_times_interval', ...
    'n',n,'q',q,'x_dimension',n+1,'ny',1,'seed',o.Seed, ...
    'eta',o.Eta,'alpha_source',o.Alpha,'orientation',orientation, ...
    'number_base_rays',number_base,'number_lifted_rays',number_base+1, ...
    'number_x_vertices',2*(number_base+1),'number_y_vertices',2, ...
    'number_pgm_active_constraints',2*n+1, ...
    'active_constraint_excess',n,'edge_degeneracy_excess',number_base-n, ...
    'old_ray_degree',q*(n-q)+1,'z_ray_degree',number_base, ...
    'householder_reflections',o.HouseholderReflections, ...
    'scaling_range',o.ScalingRange,'nonnegative_margin',o.NonnegativeMargin, ...
    'materialized_x_vertices',materialize);
instance.certificate=struct( ...
    'Sstar',Sstar,'m',m,'M',M,'D',D,'w',w, ...
    'pgm_source_value',0,'global_source_value',-1, ...
    'incumbent_source_value',o.Alpha,'objective_offset',kappa, ...
    'pgm_local_margin',1,'unique_global_minimizer',true, ...
    'lambda_z',(3-o.Alpha)/2);
instance.generation=struct( ...
    'Mx',Mx,'My',My,'shift_x',shift_x,'shift_y',shift_y, ...
    'Q0',Q0,'c0',c0,'d0',d0, ...
    'A0',A0,'rhs_A0',rhs_A0,'B0',B0,'rhs_B0',rhs_B0, ...
    'Vx0',Vx0,'Vy0',Vy0,'subset_labels',labels, ...
    'xbar0',xbar0,'ybar0',ybar0,'xglobal0',xglobal0,'yglobal0',yglobal0, ...
    'xalpha0',xalpha0,'yalpha0',yalpha0,'zalpha',zalpha);

if o.Validate
    instance.validation=validate_lifted_highdim_dblp(instance);
    if ~instance.validation.passed
        error('generate_lifted_highdim_dblp:ValidationFailed','%s',instance.validation.message);
    end
else
    instance.validation=struct('passed',NaN,'message','Validation was not requested.');
end
save_file=char(o.SaveFile);
if ~isempty(save_file), save(save_file,'instance','-v7.3'); end
end

function materialize=resolve_materialization(request,elements,limit)
if islogical(request)
    materialize=request;
elseif strcmpi(char(request),'auto')
    materialize=elements<=limit;
else
    error('generate_lifted_highdim_dblp:InvalidMaterialization', ...
        'MaterializeXVertices must be true, false, or ''auto''.');
end
end
