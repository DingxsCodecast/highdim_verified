function results = run_jgo_solver_experiments(output_directory)
%RUN_JGO_SOLVER_EXPERIMENTS Formula-blind solver diagnostics for the lift.
%   Uses the displayed H-representations and objective coefficients only.
%   No planted solution or analytic ray-depth formula enters either solver.

if nargin < 1 || isempty(output_directory)
    output_directory = fullfile(fileparts(mfilename('fullpath')), ...
        'results','jgo_202609');
end
if ~isfolder(output_directory), mkdir(output_directory); end

configs = {4,2,'none',11; 6,2,'HH',12; 6,3,'QR',13; 8,3,'QR',14};
starts = 0:0.1:1;
lp_options = optimoptions('linprog','Display','none');
nlp_options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',500,'OptimalityTolerance',1e-9, ...
    'ConstraintTolerance',1e-9,'StepTolerance',1e-12);
raw = cell(size(configs,1)*numel(starts)*2,14);
exhaustive = cell(3,11);
cursor = 0; ecursor = 0;

for k = 1:size(configs,1)
    n = configs{k,1}; q = configs{k,2}; orientation = configs{k,3};
    seed = configs{k,4};
    instance = generate_lifted_highdim_dblp(n,'Q',q,'Eta',0.5, ...
        'Alpha',-0.5,'Orientation',orientation,'Seed',seed, ...
        'MaterializeXVertices',true,'Validate',true);
    assert(instance.validation.passed);
    x0 = 0.8*instance.pgm_x + 0.2*mean(instance.Vx,2);
    for j = 1:numel(starts)
        y0_source = starts(j);
        y0 = instance.generation.My*y0_source + instance.generation.shift_y;

        tic_lp = tic;
        [x_lp,y_lp,lp_flag,lp_iterations] = alternating_lp(instance,y0,lp_options);
        lp_seconds = toc(tic_lp);
        cursor = cursor + 1;
        raw(cursor,:) = make_row(instance,orientation,seed,'alternating LP', ...
            y0_source,x_lp,y_lp,lp_flag,lp_iterations,lp_seconds);

        tic_nlp = tic;
        [z_nlp,~,nlp_flag,nlp_output] = fmincon( ...
            @(z)objective(instance,z),[x0;y0], ...
            blkdiag(instance.A,instance.B), ...
            [instance.rhs_A;instance.rhs_B],[],[],[],[],[],nlp_options);
        nlp_seconds = toc(tic_nlp);
        cursor = cursor + 1;
        raw(cursor,:) = make_row(instance,orientation,seed,'fmincon SQP', ...
            y0_source,z_nlp(1:end-1),z_nlp(end),nlp_flag, ...
            nlp_output.iterations,nlp_seconds);
    end

    if n <= 6
        tic_enum = tic;
        [VX,xinfo] = enumerate_hrep_vertices(instance.A,instance.rhs_A);
        [VY,yinfo] = enumerate_hrep_vertices(instance.B,instance.rhs_B);
        values = VX.'*instance.Q*VY + (instance.c.'*VX).'*ones(1,size(VY,2)) ...
            + ones(size(VX,2),1)*(instance.d.'*VY);
        [best,linear_index] = min(values,[],'all','linear');
        [ix,iy] = ind2sub(size(values),linear_index);
        ecursor = ecursor + 1;
        exhaustive(ecursor,:) = {n,q,string(orientation),seed,size(VX,2),size(VY,2), ...
            xinfo.candidate_bases,yinfo.candidate_bases,best-instance.minimum, ...
            max(norm(VX(:,ix)-instance.solution_x,inf), ...
                norm(VY(:,iy)-instance.solution_y,inf)),toc(tic_enum)};
    end
    fprintf('JGO solver diagnostic completed: n=%d q=%d %s seed=%d\n', ...
        n,q,orientation,seed);
end

raw = cell2table(raw(1:cursor,:), 'VariableNames', ...
    {'n','q','orientation','seed','method','y0_source','value', ...
     'objective_gap','maximum_violation','exitflag','iterations', ...
     'seconds','at_pgm','at_global'});
exhaustive = cell2table(exhaustive(1:ecursor,:), 'VariableNames', ...
    {'n','q','orientation','seed','x_vertices','y_vertices', ...
     'x_candidate_bases','y_candidate_bases','objective_error', ...
     'solution_error','seconds'});
summary = groupsummary(raw,{'n','q','orientation','method'}, ...
    {'sum','mean','max'}, {'at_pgm','at_global','objective_gap', ...
    'maximum_violation','seconds'});
assert(all(raw.maximum_violation < 1e-6));
assert(all(exhaustive.objective_error < 1e-7));
assert(all(exhaustive.solution_error < 1e-6));

writetable(raw,fullfile(output_directory,'solver_runs.csv'));
writetable(summary,fullfile(output_directory,'solver_summary.csv'));
writetable(exhaustive,fullfile(output_directory,'hrep_global_checks.csv'));
results = struct('raw',raw,'summary',summary,'exhaustive',exhaustive, ...
    'matlab_version',version,'configurations',configs,'y_starts',starts, ...
    'linprog_options',lp_options,'fmincon_options',nlp_options);
save(fullfile(output_directory,'solver_results.mat'),'results','-v7.3');
fprintf('JGO solver diagnostics saved in %s\n',output_directory);
end

function [x,y,exitflag,iterations] = alternating_lp(instance,y,options)
exitflag = 1;
iterations = 0;
for k = 1:20
    [x,~,flag_x] = linprog(instance.Q*y+instance.c, ...
        instance.A,instance.rhs_A,[],[],[],[],options);
    if flag_x <= 0
        error('run_jgo_solver_experiments:LPFailure','x-block LP failed.');
    end
    [y_new,~,flag_y] = linprog(instance.Q.'*x+instance.d, ...
        instance.B,instance.rhs_B,[],[],[],[],options);
    if flag_y <= 0
        error('run_jgo_solver_experiments:LPFailure','y-block LP failed.');
    end
    iterations = k;
    if abs(y_new-y) <= 1e-10*(1+abs(y))
        y = y_new;
        return;
    end
    y = y_new;
end
exitflag = 0;
end

function f = objective(instance,z)
x = z(1:end-1); y = z(end);
f = x.'*instance.Q*y + instance.c.'*x + instance.d.'*y;
end

function row = make_row(instance,orientation,seed,method,y0_source, ...
    x,y,exitflag,iterations,seconds)
value = x.'*instance.Q*y + instance.c.'*x + instance.d.'*y;
gap = value-instance.minimum;
violation = max([instance.A*x-instance.rhs_A; ...
    instance.B*y-instance.rhs_B;0]);
at_pgm = norm(x-instance.pgm_x,inf) < 1e-5 && ...
    abs(y-instance.pgm_y) < 1e-5;
at_global = norm(x-instance.solution_x,inf) < 1e-5 && ...
    abs(y-instance.solution_y) < 1e-5;
row = {instance.metadata.n,instance.metadata.q,string(upper(orientation)), ...
    seed,string(method),y0_source,value,gap,violation,exitflag, ...
    iterations,seconds,at_pgm,at_global};
end
