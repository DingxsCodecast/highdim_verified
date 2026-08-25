function results=stress_test_highdim_generator()
%STRESS_TEST_HIGHDIM_GENERATOR Cross dimensions, q, seeds, and orientations.
addpath(fileparts(mfilename('fullpath')));
configs=[4 2;5 2;8 2;12 2;25 2;50 2;100 2;200 2;500 2;6 3;8 3;12 4];
seeds=[1 7 2026]; orientations={'HH','QR'};
records=cell(size(configs,1)*numel(seeds)*numel(orientations)*2,8); cursor=0;
for r=1:size(configs,1)
    n=configs(r,1); q=configs(r,2);
    for seed=seeds
        for oi=1:numel(orientations)
            orientation=orientations{oi};
            for mode=1:2
                cursor=cursor+1; timer=tic;
                if mode==1
                    inst=generate_highdim_dblp(n,'Q',q,'Seed',seed,'Orientation',orientation, ...
                        'YDimension',max(1,ceil(n/2)),'MaterializeXVertices',n<=50,'Validate',true);
                    mode_name='structural';
                else
                    inst=generate_lifted_highdim_dblp(n,'Q',q,'Seed',seed,'Orientation',orientation, ...
                        'Eta',0.25+0.5*mod(seed,3),'Alpha',-0.25-0.25*mod(seed,3), ...
                        'MaterializeXVertices',n<=50,'Validate',true);
                    mode_name='lifted';
                end
                records(cursor,:)={string(mode_name),n,q,seed,string(orientation),toc(timer),inst.validation.passed,nchoosek(n,q)};
            end
        end
    end
end
table_records=cell2table(records,'VariableNames',{'mode','n','q','seed','orientation','seconds','passed','rays'});
results=struct('passed',all(table_records.passed),'cases',height(table_records), ...
    'records',table_records,'maximum_seconds',max(table_records.seconds));
assert(results.passed);
fprintf('Stress test passed: %d cases; maximum %.3f s.\n',results.cases,results.maximum_seconds);
end
