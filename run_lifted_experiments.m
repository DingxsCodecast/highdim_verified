function results = run_lifted_experiments(output_directory, varargin)
%RUN_LIFTED_EXPERIMENTS Reproduce the paper's structural and cut experiments.

if nargin<1||isempty(output_directory)
    output_directory=fullfile(fileparts(mfilename('fullpath')),'results');
end
parser=inputParser;
addParameter(parser,'Profile','full',@(x)ischar(x)||isstring(x));
parse(parser,varargin{:});
profile=lower(char(parser.Results.Profile));
if ~ismember(profile,{'smoke','full'})
    error('run_lifted_experiments:InvalidProfile','Profile must be smoke or full.');
end
if ~isfolder(output_directory), mkdir(output_directory); end

tol_feas=1e-9; tol_cut=1e-9; tol_cone=1e-9;
fprintf('Phase 0: analytic smoke test...\n');
smoke=generate_lifted_highdim_dblp(4,'Q',2,'Eta',0.5,'Alpha',-0.5, ...
    'Orientation','none','Seed',20260820,'Validate',true);
tri_smoke=reference_triangulation(4,2,'CoverageSamples',200,'Seed',20260820);
assert(tri_smoke.passed,tri_smoke.message);
star_cell=find(any(tri_smoke.cells==find(ismember(tri_smoke.labels,smoke.certificate.Sstar,'rows')),2),1);
smoke_cut=oracle_polar_cut(smoke,tri_smoke,star_cell,'Tolerance',tol_cut);
assert(smoke_cut.passed);
smoke_depths=generic_ray_depths(smoke,tri_smoke.labels);
smoke_numerical_cut=numerical_polar_cut(smoke,tri_smoke,star_cell,smoke_depths, ...
    'Tolerance',tol_cut);
assert(smoke_depths.passed&&smoke_numerical_cut.passed);
assert(abs(smoke.pgm_value)<=tol_feas&&abs(smoke.incumbent_value+0.5)<=tol_feas&&abs(smoke.minimum+1)<=tol_feas);
assert(abs(smoke.certificate.lambda_z-1.75)<=tol_cut);
assert(abs(smoke_cut.incumbent_lhs-1)<=tol_cut);
assert(abs(smoke_cut.global_lhs-4/3.5)<=tol_cut);

fprintf('Phase 1-2: structural and lifted neighborhood diagnostics...\n');
if strcmp(profile,'smoke')
    neighborhood_configs=[4 2;8 3;10 5];
else
    neighborhood_configs=[4 2;5 2;8 2;10 2;20 2;50 2;100 2;200 2;500 2; ...
        6 3;8 3;10 3;15 3;20 3;30 3;50 3;8 4;10 5;12 6];
end
neighbor_rows=cell(size(neighborhood_configs,1),13);
for r=1:size(neighborhood_configs,1)
    n=neighborhood_configs(r,1); q=neighborhood_configs(r,2);
    labels=deterministic_subset_sample(n,q,min(50,nchoosek(n,q)));
    exact_times=zeros(size(labels,1),1); baseline_times=exact_times;
    exact=true(size(labels,1),1); baseline_exact=exact;
    for k=1:size(labels,1)
        timer=tic; truth=delta_neighbors(labels(k,:),n); exact_times(k)=toc(timer);
        timer=tic; returned=baseline_nminus1_neighbors(labels(k,:),n); baseline_times(k)=toc(timer);
        exact(k)=size(truth,1)==q*(n-q);
        baseline_exact(k)=isequal(returned,truth);
    end
    degree=q*(n-q); deficit=degree-(n-1);
    neighbor_rows(r,:)={n,q,size(labels,1),degree,n-1,deficit,(n-1)/degree, ...
        all(exact),all(baseline_exact),mean(exact_times),mean(baseline_times),degree+1,nchoosek(n,q)};
end
neighborhood=cell2table(neighbor_rows,'VariableNames', ...
    {'n','q','sampled_rays','true_degree','simple_degree','deficit','baseline_recall', ...
     'oracle_exact','baseline_exact','oracle_mean_seconds','baseline_mean_seconds', ...
     'lifted_old_ray_degree','lifted_z_ray_degree'});

fprintf('Phase 2b: independent H-description vertex and adjacency recovery...\n');
if strcmp(profile,'smoke')
    hrep_configs=[4 2;6 3];
else
    hrep_configs=[4 2;5 2;6 2;7 2;8 2;6 3;7 3;8 3];
end
hrep_rows=cell(size(hrep_configs,1),16);
for r=1:size(hrep_configs,1)
    n=hrep_configs(r,1); q=hrep_configs(r,2);
    [AH,bH]=hypersimplex_hrep(n,q);
    timer=tic; [VH,enum_info]=enumerate_hrep_vertices(AH,bH,'Tolerance',tol_feas);
    enumeration_seconds=toc(timer);
    timer=tic; graph=hrep_adjacency_graph(AH,bH,VH,'Tolerance',tol_feas);
    adjacency_seconds=toc(timer);
    labels=enumerate_q_subsets(n,q);
    H=zeros(size(labels,1),n-1);
    for k=1:size(labels,1), H(k,:)=subset_to_h(labels(k,:),n).'; end
    [H,order]=sortrows(H,1:n-1); labels=labels(order,:);
    vertex_error=max(abs(VH.'-H),[],'all');
    exchange=false(size(labels,1));
    for k=1:size(labels,1)
        neighbors=delta_neighbors(labels(k,:),n);
        [found,loc]=ismember(neighbors,labels,'rows');
        assert(all(found)); exchange(k,loc)=true;
    end
    expected_edges=nchoosek(n,q)*q*(n-q)/2;
    hrep_rows(r,:)={n,q,size(AH,1),enum_info.candidate_bases,size(VH,2),nchoosek(n,q), ...
        vertex_error<=100*tol_feas,vertex_error,nnz(triu(graph.adjacency)),expected_edges, ...
        isequal(graph.adjacency,exchange),min(graph.degree),max(graph.degree),q*(n-q), ...
        enumeration_seconds,adjacency_seconds};
end
independent_neighborhood=cell2table(hrep_rows,'VariableNames', ...
    {'n','q','inequalities','candidate_bases','hrep_vertices','expected_vertices', ...
     'vertex_sets_match','max_vertex_error','hrep_edges','expected_edges', ...
     'edge_sets_match','minimum_degree','maximum_degree','expected_degree', ...
     'enumeration_seconds','adjacency_seconds'});
witness=porembski_first_step_witness('Tolerance',tol_feas);
procedure_witness=table(witness.dimension,witness.enumerated_vertices, ...
    witness.actual_neighbors,witness.n_minus_one_neighbors,witness.neighbor_excess, ...
    witness.possible_n_minus_one_selections,witness.first_step_count_matches,witness.passed, ...
    'VariableNames',{'dimension','enumerated_vertices','actual_neighbors', ...
    'n_minus_one_neighbors','neighbor_excess','possible_n_minus_one_selections', ...
    'first_step_count_matches','witness_passed'});

fprintf('Phase 3: reference triangulation and coverage...\n');
if strcmp(profile,'smoke'), tri_configs=[4 2;6 3];
else, tri_configs=[4 2;5 2;6 2;7 2;8 2;6 3;7 3;8 3]; end
triangulations=cell(size(tri_configs,1),1);
decomp_rows=cell(size(tri_configs,1),11);
for r=1:size(tri_configs,1)
    n=tri_configs(r,1); q=tri_configs(r,2);
    timer=tic; tri=reference_triangulation(n,q,'CoverageSamples',500,'Seed',1000+n*10+q); total=toc(timer);
    assert(tri.passed,tri.message); triangulations{r}=tri;
    decomp_rows(r,:)={n,q,size(tri.H,1),tri.number_cells,all(tri.sigma_min>tri.rank_tolerances), ...
        tri.coverage_rate,tri.volume_sum,tri.analytic_volume,tri.volume_relative_error, ...
        tri.discarded_rank_deficient,total};
end
decomposition=cell2table(decomp_rows,'VariableNames', ...
    {'n','q','input_rays','cells','all_simplicial','coverage_rate','volume_sum', ...
     'analytic_volume','volume_relative_error','discarded_rank_deficient','seconds'});

fprintf('Phase 4: black-box ray extension and numerical polar-cut validity...\n');
if strcmp(profile,'smoke'), etas=0.5; alphas=-0.5;
else, etas=[0.25 0.5 1.0 1.5]; alphas=[-0.25 -0.5 -0.75]; end
cut_rows=cell(numel(triangulations)*numel(etas)*numel(alphas),22); cursor=0;
for r=1:numel(triangulations)
    tri=triangulations{r}; n=tri.n; q=tri.q;
    for eta=etas
        for alpha=alphas
            cursor=cursor+1;
            instance=generate_lifted_highdim_dblp(n,'Q',q,'Eta',eta,'Alpha',alpha, ...
                'Orientation','none','Seed',2000+n*10+q,'Validate',true);
            depths=generic_ray_depths(instance,tri.labels,'Coordinates','source');
            exact_lambda=zeros(size(depths.lambda));
            for j=1:size(tri.labels,1)
                exact_lambda(j)=lambda_oracle(instance,tri.labels(j,:));
            end
            exact_lambda(end)=lambda_oracle(instance,'Z');
            depth_error=abs(depths.lambda-exact_lambda)./exact_lambda;
            assert(depths.passed&&max(depth_error)<=100*tol_cut);
            timer=tic; max_cut_error=0; global_count=0; global_ok=true; incumbent_ok=true;
            dmins=zeros(tri.number_cells,1); dmeans=dmins; dgeos=dmins; logs=dmins;
            for k=1:tri.number_cells
                cut=numerical_polar_cut(instance,tri,k,depths,'Tolerance',tol_cut);
                assert(cut.passed);
                max_cut_error=max(max_cut_error,cut.max_cut_recovery_error);
                dmins(k)=cut.D_min; dmeans(k)=cut.D_mean; dgeos(k)=cut.D_geo; logs(k)=cut.log_volume_multiplier;
                if cut.contains_global
                    global_count=global_count+1;
                    global_ok=global_ok&&(cut.global_lhs>=1-tol_cut)&& ...
                        abs(cut.global_lhs-4/(3-alpha))<=100*tol_cut;
                end
                if cut.contains_incumbent, incumbent_ok=incumbent_ok&&abs(cut.incumbent_lhs-1)<=100*tol_cut; end
            end
            cut_seconds=toc(timer);
            cut_rows(cursor,:)={n,q,eta,alpha,numel(depths.lambda),tri.number_cells, ...
                max(depth_error),mean(depth_error),max(abs(depths.residual)), ...
                sum(depths.function_evaluations),max_cut_error,global_count,global_ok, ...
                incumbent_ok,min(dmins),mean(dmeans),exp(mean(log(dgeos))),mean(logs), ...
                depths.seconds,cut_seconds,tol_cut,tol_cone};
        end
    end
end
cut_depth=cell2table(cut_rows,'VariableNames', ...
    {'n','q','eta','alpha','rays','cells','max_blackbox_relative_depth_error', ...
     'mean_blackbox_relative_depth_error','max_blackbox_level_residual', ...
     'blackbox_function_evaluations','max_cut_recovery_error','global_cells', ...
     'global_preserved','incumbent_on_boundary','D_min_min','D_mean_mean', ...
     'D_geo_geomean','mean_log_volume_multiplier','blackbox_seconds', ...
     'cut_assembly_seconds','tol_cut','tol_cone'});

fprintf('Phase 5: none/HH/QR orientation robustness...\n');
if strcmp(profile,'smoke'), orient_configs=[4 2]; seeds=1;
else, orient_configs=[6 2;6 3;8 3]; seeds=1:3; end
orientations={'none','HH','QR'};
orientation_rows=cell(size(orient_configs,1)*numel(seeds)*numel(orientations),15); cursor=0;
for r=1:size(orient_configs,1)
    n=orient_configs(r,1); q=orient_configs(r,2);
    tri=reference_triangulation(n,q,'CoverageSamples',200,'Seed',3000+n*10+q);
    for seed=seeds
        for oidx=1:numel(orientations)
            orientation=orientations{oidx}; cursor=cursor+1;
            timer=tic;
            instance=generate_lifted_highdim_dblp(n,'Q',q,'Eta',0.5,'Alpha',-0.5, ...
                'Orientation',orientation,'Seed',seed,'Validate',true);
            generation_seconds=toc(timer);
            depths=generic_ray_depths(instance,tri.labels,'Coordinates','transformed');
            exact_lambda=zeros(size(depths.lambda));
            for j=1:size(tri.labels,1)
                exact_lambda(j)=lambda_oracle(instance,tri.labels(j,:));
            end
            exact_lambda(end)=lambda_oracle(instance,'Z');
            depth_error=abs(depths.lambda-exact_lambda)./exact_lambda;
            max_cut_error=0; global_ok=true; cut_timer=tic;
            for k=1:tri.number_cells
                cut=numerical_polar_cut(instance,tri,k,depths, ...
                    'Coordinates','transformed','Tolerance',tol_cut);
                max_cut_error=max(max_cut_error,cut.max_cut_recovery_error);
                if cut.contains_global, global_ok=global_ok&&cut.global_lhs>=1-tol_cut; end
            end
            cut_seconds=toc(cut_timer);
            orientation_rows(cursor,:)={n,q,seed,string(upper(orientation)),tri.number_cells,tri.passed, ...
                instance.validation.passed,max(depth_error),max(abs(depths.residual)), ...
                max_cut_error,global_ok,generation_seconds,depths.seconds,cut_seconds, ...
                instance.validation.metrics.maximum_transformed_boundary_residual};
        end
    end
end
orientation_robustness=cell2table(orientation_rows,'VariableNames', ...
    {'n','q','seed','orientation','cells','decomposition_passed','instance_passed', ...
     'max_blackbox_relative_depth_error','max_blackbox_level_residual', ...
     'max_cut_recovery_error','global_preserved','generation_seconds', ...
     'blackbox_seconds','cut_assembly_seconds','max_transformed_boundary_residual'});

fprintf('Generation timings and deterministic search-to-cut pipeline...\n');
if strcmp(profile,'smoke')
    generation_configs={'structural',4,2,'HH';'lifted',4,2,'HH'}; timing_seeds=1:2;
else
    generation_configs={ ...
        'structural',4,2,'HH';'structural',20,2,'HH';'structural',100,2,'HH'; ...
        'structural',500,2,'HH';'structural',20,3,'HH';'structural',20,3,'QR'; ...
        'lifted',4,2,'HH';'lifted',20,2,'HH';'lifted',100,2,'HH'; ...
        'lifted',500,2,'HH';'lifted',20,3,'HH';'lifted',20,3,'QR'};
    timing_seeds=1:5;
end
generation_rows=cell(size(generation_configs,1),8);
for r=1:size(generation_configs,1)
    mode=generation_configs{r,1}; n=generation_configs{r,2}; q=generation_configs{r,3}; orientation=generation_configs{r,4};
    times=zeros(numel(timing_seeds),1); passes=false(size(times));
    for k=1:numel(timing_seeds)
        timer=tic;
        if strcmp(mode,'structural')
            inst=generate_highdim_dblp(n,'Q',q,'Orientation',orientation,'Seed',timing_seeds(k),'Validate',true);
        else
            inst=generate_lifted_highdim_dblp(n,'Q',q,'Orientation',orientation,'Seed',timing_seeds(k),'Validate',true);
        end
        times(k)=toc(timer); passes(k)=inst.validation.passed;
    end
    generation_rows(r,:)={string(mode),n,q,string(orientation),numel(times),median(times),max(times),all(passes)};
end
generation=cell2table(generation_rows,'VariableNames', ...
    {'mode','n','q','orientation','seeds','median_seconds','maximum_seconds','all_passed'});

pipeline_configs=[4 2;6 2;6 3];
pipeline_rows=cell(size(pipeline_configs,1),10);
for r=1:size(pipeline_configs,1)
    n=pipeline_configs(r,1); q=pipeline_configs(r,2);
    timer=tic;
    inst=generate_lifted_highdim_dblp(n,'Q',q,'Orientation','QR','Seed',4000+r,'Validate',true);
    tri=reference_triangulation(n,q,'CoverageSamples',200,'Seed',4000+r);
    depths=generic_ray_depths(inst,tri.labels,'Coordinates','transformed');
    all_cuts=true; global_ok=true;
    for k=1:tri.number_cells
        cut=numerical_polar_cut(inst,tri,k,depths,'Coordinates','transformed','Tolerance',tol_cut);
        all_cuts=all_cuts&&cut.passed;
        if cut.contains_global, global_ok=global_ok&&cut.global_lhs>=1-tol_cut; end
    end
    pipeline_rows(r,:)={n,q,0,true,true,tri.number_cells,all_cuts,global_ok,inst.minimum,toc(timer)};
end
pipeline=cell2table(pipeline_rows,'VariableNames', ...
    {'n','q','initial_y_source','PGM_reached','PGM_matches_planted','child_cells', ...
     'all_cuts_valid','global_preserved','known_global_value_transformed','total_seconds'});

environment=capture_environment();
results=struct('profile',profile,'created',char(datetime('now','TimeZone','local')), ...
    'tolerances',struct('feasibility',tol_feas,'cut',tol_cut,'cone',tol_cone), ...
    'smoke',smoke.validation,'smoke_cut',smoke_cut,'smoke_depths',smoke_depths, ...
    'smoke_numerical_cut',smoke_numerical_cut,'neighborhood',neighborhood, ...
    'independent_neighborhood',independent_neighborhood,'procedure_witness',procedure_witness, ...
    'decomposition',decomposition,'cut_depth',cut_depth, ...
    'orientation_robustness',orientation_robustness,'generation',generation, ...
    'pipeline',pipeline,'environment',environment);

writetable(neighborhood,fullfile(output_directory,'table_neighborhood.csv'));
writetable(independent_neighborhood,fullfile(output_directory,'table_independent_neighborhood.csv'));
writetable(procedure_witness,fullfile(output_directory,'table_procedure_witness.csv'));
writetable(decomposition,fullfile(output_directory,'table_decomposition.csv'));
writetable(cut_depth,fullfile(output_directory,'table_cut_depth.csv'));
writetable(orientation_robustness,fullfile(output_directory,'table_orientation.csv'));
writetable(generation,fullfile(output_directory,'table_generation.csv'));
writetable(pipeline,fullfile(output_directory,'table_pipeline.csv'));
save(fullfile(output_directory,'experiment_results.mat'),'results','-v7.3');
write_environment(fullfile(output_directory,'environment.txt'),environment);
fprintf('All %s-profile experiments passed. Results: %s\n',profile,output_directory);
end

function environment=capture_environment()
environment=struct();
environment.matlab_version=version;
environment.matlab_release=version('-release');
environment.computer=computer;
environment.matlab_numcores=feature('numcores');
environment.parallel_pool='none';
try
    pool=gcp('nocreate'); if ~isempty(pool), environment.parallel_pool=sprintf('%d workers',pool.NumWorkers); end
catch
    environment.parallel_pool='Parallel Computing Toolbox unavailable';
end
try
    environment.gpu_count=gpuDeviceCount('available');
catch
    environment.gpu_count=0;
end
environment.cpu=strtrim(getenv('PROCESSOR_IDENTIFIER'));
environment.logical_processors=strtrim(getenv('NUMBER_OF_PROCESSORS'));
environment.windows=strtrim(getenv('OS'));
try
    [~,system_memory]=memory;
    environment.ram_bytes=system_memory.PhysicalMemory.Total;
catch
    environment.ram_bytes=NaN;
end
end

function write_environment(file_name,environment)
fid=fopen(file_name,'w'); cleanup=onCleanup(@()fclose(fid));
names=fieldnames(environment);
for k=1:numel(names), fprintf(fid,'%s: %s\n',names{k},string(environment.(names{k}))); end
end
