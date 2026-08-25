function results=run_generator_tests()
%RUN_GENERATOR_TESTS Regression tests for structural and lifted modes.
addpath(fileparts(mfilename('fullpath')));
fprintf('1/15 structural q=2 exhaustive...\n');
s=generate_highdim_dblp(4,'Seed',11,'MaterializeXVertices',true);
assert(s.validation.passed&&s.metadata.number_apex_neighbors==6);

fprintf('2/15 structural general q and QR...\n');
s3=generate_highdim_dblp(8,'Q',3,'YDimension',5,'Orientation','QR','Seed',12);
assert(s3.validation.passed&&s3.metadata.cross_section_vertex_degree==15);

fprintf('3/15 structural implicit large mode...\n');
sl=generate_highdim_dblp(200,'Q',2,'MaterializeXVertices',false,'Seed',13);
assert(sl.validation.passed&&isempty(sl.Vx));

fprintf('4/15 lifted hand-check smoke identities...\n');
L=generate_lifted_highdim_dblp(4,'Q',2,'Eta',0.5,'Alpha',-0.5,'Orientation','none','Seed',14);
assert(L.validation.passed);
assert(max(abs([L.pgm_value,L.incumbent_value,L.minimum,L.certificate.lambda_z]-[0,-0.5,-1,1.75]))<1e-10);

fprintf('5/15 lifted general q and QR affine invariance...\n');
L3=generate_lifted_highdim_dblp(8,'Q',3,'Eta',1.5,'Alpha',-0.75,'Orientation','QR','Seed',15);
assert(L3.validation.passed&&L3.metadata.old_ray_degree==16);

fprintf('6/15 exact neighborhood oracle...\n');
S=[1 3 6]; N=delta_neighbors(S,8); B=baseline_nminus1_neighbors(S,8);
assert(size(N,1)==15&&size(unique(N,'rows'),1)==15&&size(B,1)==7&&~isequal(N,B));

fprintf('7/15 independent H-description vertex and adjacency recovery...\n');
[AH,bH]=hypersimplex_hrep(6,3);
[VH,IH]=enumerate_hrep_vertices(AH,bH);
GH=hrep_adjacency_graph(AH,bH,VH);
assert(IH.number_vertices==20&&all(GH.degree==9)&&nnz(triu(GH.adjacency))==90);

fprintf('8/15 procedure-level four-dimensional neighbor-count witness...\n');
W=porembski_first_step_witness();
assert(W.passed&&W.actual_neighbors==4&&W.n_minus_one_neighbors==3&& ...
    W.possible_n_minus_one_selections==4);

fprintf('9/15 reference triangulation and volume...\n');
T=reference_triangulation(6,3,'CoverageSamples',300,'Seed',16);
assert(T.passed&&T.number_cells>0&&T.volume_relative_error<1e-10);

fprintf('10/15 oracle cell-specific cut in source coordinates...\n');
Ls=generate_lifted_highdim_dblp(6,'Q',3,'Orientation','none','Seed',16);
C=oracle_polar_cut(Ls,T,1);
assert(C.passed&&C.max_relative_depth_error<1e-10);

fprintf('11/15 black-box ray extension and numerical cut...\n');
Ds=generic_ray_depths(Ls,T.labels);
truth=zeros(size(Ds.lambda));
for k=1:size(T.labels,1), truth(k)=lambda_oracle(Ls,T.labels(k,:)); end
truth(end)=lambda_oracle(Ls,'Z');
assert(Ds.passed&&max(abs(Ds.lambda-truth)./truth)<1e-10);
Cs=numerical_polar_cut(Ls,T,1,Ds);
assert(Cs.passed&&Cs.max_cut_recovery_error<1e-10);

fprintf('12/15 black-box cut after HH transformation...\n');
Lh=generate_lifted_highdim_dblp(6,'Q',3,'Orientation','HH','Seed',16);
Dt=generic_ray_depths(Lh,T.labels,'Coordinates','transformed');
Ct=numerical_polar_cut(Lh,T,1,Dt,'Coordinates','transformed');
assert(Dt.passed&&Ct.passed&&Ct.max_cut_recovery_error<1e-10);

fprintf('13/15 reproducibility and caller RNG preservation...\n');
rng(9876,'twister'); before=rng;
a=generate_lifted_highdim_dblp(7,'Q',3,'Seed',17,'Validate',false);
after=rng; b=generate_lifted_highdim_dblp(7,'Q',3,'Seed',17,'Validate',false);
assert(isequal(before,after)&&isequal(a.Q,b.Q)&&isequal(a.certificate.Sstar,b.certificate.Sstar));

fprintf('14/15 positional wrapper and save/load...\n');
[Q,A,B,rhs_A,rhs_B,c,d,fstar,xstar,ystar,certificate,wrapped]=highdim_generator(8,'Q',3,'Seed',18);
assert(wrapped.validation.passed&&certificate.unique_global_minimizer);
assert(abs(xstar'*Q*ystar+c'*xstar+d'*ystar-fstar)<1e-8);
assert(size(A,1)==numel(rhs_A)&&size(B,1)==numel(rhs_B));

fprintf('15/15 invalid options fail clearly...\n');
assert_throws(@()generate_highdim_dblp(4,'Q',3));
assert_throws(@()generate_lifted_highdim_dblp(8,'Eta',2));
assert_throws(@()generate_lifted_highdim_dblp(8,'Alpha',0));
assert_throws(@()generate_highdim_dblp(8,'Orientation','mystery'));

results=struct('passed',true,'tests',15,'message','All 15 generator regression tests passed.');
fprintf('%s\n',results.message);
end

function assert_throws(action)
did_throw=false;
try
    action();
catch
    did_throw=true;
end
assert(did_throw);
end
