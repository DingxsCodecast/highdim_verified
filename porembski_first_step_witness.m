function witness = porembski_first_step_witness(varargin)
%POREMBSKI_FIRST_STEP_WITNESS Procedure-level count failure on Delta(4,2).
%   The ray labelled {1,2} is recovered from the H-description and has four
%   neighbors, while an n-1 rule requests three and leaves four possible
%   three-neighbor selections.

[A,rhs] = hypersimplex_hrep(4,2);
[vertices,enumeration] = enumerate_hrep_vertices(A,rhs,varargin{:});
graph = hrep_adjacency_graph(A,rhs,vertices,varargin{:});
target = subset_to_h([1 2],4);
[distance,index] = min(max(abs(vertices-target),[],1));
if distance > 1e-8
    error('porembski_first_step_witness:TargetNotFound', ...
        'The independently enumerated H-vertices do not contain ray {1,2}.');
end
actual = graph.degree(index);
required = 3;
witness = struct('dimension',4,'ray_label',[1 2], ...
    'enumerated_vertices',size(vertices,2),'actual_neighbors',actual, ...
    'n_minus_one_neighbors',required,'neighbor_excess',actual-required, ...
    'possible_n_minus_one_selections',nchoosek(actual,required), ...
    'first_step_count_matches',actual==required,'enumeration',enumeration, ...
    'passed',size(vertices,2)==6&&actual==4&&required==3);
end
