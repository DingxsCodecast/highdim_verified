# Verifiable High-Dimensional Degenerate DBLP Benchmarks

This archive accompanies the article **"Certified High-Dimensional Testbeds for Degenerate Disjoint Bilinear Optimization"**, prepared for submission to *Journal of Global Optimization*.

Authors: Xi Chen, Xiuming Li, Xiangqi Tai, and Xiaosong Ding  
Corresponding author: Xiaosong Ding, International Business School, Beijing Foreign Studies University, Beijing, China  
Email: xiaosong.ding@hotmail.com

The code provides two benchmark modes built on the same hypersimplex geometry, `Delta(n,q)`:

- **Structural mode:** `X_(0,q)` is a pyramid over `Delta(n,q)`. Its degenerate apex is the unique global optimizer. This mode is intended for neighborhood and geometric verification.
- **Lifted mode:** `Xtilde_(0,q) = X_(0,q) x [0,1]`. Its lifted apex is a strict non-global pseudo-global minimizer (PGM), while a different vertex is the known unique global optimizer. The instance also provides a feasible incumbent with prescribed value `alpha` and closed-form maximal multipliers for all incident rays. This mode is intended for decomposition and polar-cut experiments.

Valid parameters are `n >= 4` and `2 <= q <= n-2`. Lifted instances additionally require `0 < eta < 2` and `-1 < alpha < 0`.

## 1. Extract Online Resource 1

Extract `Online_Resource_1_highdim_verified.zip` so that the
`highdim_verified` directory and this README remain together. For a MATLAB
script stored in that directory, use the directory-relative setup below:

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);
```

No machine-specific absolute path is required.

## 2. Quick start

Run both benchmark modes and a reference cut example:

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);
example_generate
```

### Structural mode

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);

S = generate_highdim_dblp(20, ...
    'Q', 3, ...
    'YDimension', 20, ...
    'Seed', 20260820, ...
    'BilinearStrength', 0.5, ...
    'Orientation', 'QR', ... % 'none', 'HH', or 'QR'
    'Validate', true);

disp(S.validation.message)
disp(S.metadata)
```

The compatibility wrapper for the original solver interface is:

```matlab
[Q,A,B,rhs_A,rhs_B,c,d,fstar,xstar,ystar,certificate,instance] = ...
    highdim_generator(20,'Q',3,'Seed',20260820);
```

### Lifted mode

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);

L = generate_lifted_highdim_dblp(8, ...
    'Q', 3, ...
    'Eta', 0.5, ...
    'Alpha', -0.5, ...
    'Orientation', 'HH', ...
    'Seed', 20260820, ...
    'Validate', true);

fprintf('PGM value       = %.16g\n',L.pgm_value);
fprintf('Incumbent value = %.16g\n',L.incumbent_value);
fprintf('Global value    = %.16g\n',L.minimum);
disp(L.certificate.Sstar)
```

With `Orientation='none'`, no scaling, rotation, or translation is applied, so the source objective values are exactly `0`, `alpha`, and `-1`. In the `HH` and `QR` modes, the objective is shifted by the constant stored in `L.certificate.objective_offset`; objective differences and canonical ray multipliers are unchanged.

Canonical rays use the vector from the PGM to an adjacent vertex, without unit-length normalization:

```matlab
Slabel = L.certificate.Sstar;
[r0,rT] = canonical_lifted_ray(L,Slabel);
lambdaS = lambda_oracle(L,Slabel);
lambdaZ = lambda_oracle(L,'Z');
```

## 3. Independent neighborhood, triangulation, and cell-cut checks

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);

N = delta_neighbors([1 3 5],8);                 % exact exchange-one-index oracle
B = baseline_nminus1_neighbors([1 3 5],8);     % intentionally incomplete diagnostic, not CDP

% Recover vertices and adjacencies from the H-description without using
% subset labels or the exchange rule.
[AH,bH] = hypersimplex_hrep(8,3);
[VH,enumInfo] = enumerate_hrep_vertices(AH,bH);
GH = hrep_adjacency_graph(AH,bH,VH);
assert(enumInfo.number_vertices == nchoosek(8,3))
assert(all(GH.degree == 3*(8-3)))

% Procedure-level neighbor-count witness on Delta(4,2).
W = porembski_first_step_witness();
assert(W.passed && W.actual_neighbors == 4 && W.n_minus_one_neighbors == 3)

T = reference_triangulation(8,3, ...
    'CoverageSamples',500,'Seed',20260820);
assert(T.passed,T.message)

L0 = generate_lifted_highdim_dblp(8,'Q',3, ...
    'Eta',0.5,'Alpha',-0.5,'Orientation','none');
C = oracle_polar_cut(L0,T,1,'Coordinates','source');
assert(C.passed)

% Black-box bracket-and-bisection depths; numerical_polar_cut does not call
% the closed-form lambda_oracle.
D = generic_ray_depths(L0,T.labels,'Coordinates','source');
CN = numerical_polar_cut(L0,T,1,D,'Coordinates','source');
assert(D.passed && CN.passed)
```

`reference_triangulation` is intended only for small and medium instances. It materializes all `nchoosek(n,q)` hypersimplex rays, asks Qhull with option `QJ` for a candidate triangulation, and independently checks cell ranks, Eulerian-number volume, and sampled-point coverage. Large instances should use certified ray sampling or an implicit ray/column-generation method.

A cell cut is branch-specific. The inequality `a'*(x-xbar) >= 1` must be combined with the current simplicial-cell constraint; a cut derived for one cell must not be imposed indiscriminately on other cells.

## 4. Regression and stress tests

In an open MATLAB graphical session, add the extracted folder to the path and run:

```matlab
addpath('PATH_TO_EXTRACTED_HIGH_DIM_VERIFIED_FOLDER')
run_generator_tests
stress_test_highdim_generator
```

The regression suite contains 15 checks covering general `q`, HH/QR transformations, hand-verifiable identities, the exchange oracle, independent H-description vertex and adjacency recovery, the four-dimensional procedure witness, black-box ray extension, analytic volume coverage, source/transformed numerical cuts, RNG reproducibility, the compatibility wrapper, and invalid inputs.

The stress suite crosses 12 `(n,q)` settings, 3 seeds, 2 orientations, and 2 modes, for 144 instances.

## 5. Reproduce the reported experiments

```matlab
repo = fileparts(mfilename('fullpath'));
addpath(repo);

out = fullfile(repo,'results');
run_lifted_experiments(out,'Profile','full');
run_jgo_solver_experiments(fullfile(repo,'results','jgo_202609'));
```

The archived outputs include:

- `experiment_results.mat`: complete MATLAB structures and raw metrics;
- `table_neighborhood.csv`: exact oracle versus the `n-1` diagnostic;
- `table_independent_neighborhood.csv`: vertices and edges independently recovered from the H-description;
- `table_procedure_witness.csv`: the actual neighbor-count witness on `Delta(4,2)`;
- `table_decomposition.csv`: cell counts, ranks, analytic volume, and coverage;
- `table_cut_depth.csv`: black-box bisection depths, analytic errors, and numerical-cut validity for 96 configurations;
- `table_orientation.csv`: affine robustness for none/HH/QR orientations and three seeds;
- `table_generation.csv`: generation and validation time over five seeds for both modes;
- `table_pipeline.csv`: deterministic search-to-cut pipeline;
- `environment.txt`: CPU, memory, operating system, MATLAB version, parallel-pool state, and GPU information.
- `results/jgo_202609/solver_runs.csv`: 88 formula-blind alternating-LP and SQP runs;
- `results/jgo_202609/solver_summary.csv`: solver outcomes by geometry and orientation;
- `results/jgo_202609/hrep_global_checks.csv`: small-scale global solutions recovered only from the displayed H-descriptions.

The reported experiment tolerances are `tol_feas = tol_cut = tol_cone = 1e-9`. Rank tests use `max(size(R))*eps(norm(R,2))`. Failure records retain subset or cell labels, cut coefficients, and the relevant residuals rather than only a status flag.

## 6. Main files

- `generate_highdim_dblp.m`, `validate_highdim_dblp.m`: structural mode;
- `generate_lifted_highdim_dblp.m`, `validate_lifted_highdim_dblp.m`: lifted mode;
- `delta_neighbors.m`: exact neighborhood oracle;
- `hypersimplex_hrep.m`, `enumerate_hrep_vertices.m`, `hrep_adjacency_graph.m`: independent H-description verification;
- `porembski_first_step_witness.m`: four-dimensional procedure-level neighbor-count witness;
- `reference_triangulation.m`: small-instance reference decomposition;
- `lambda_oracle.m`, `oracle_polar_cut.m`: analytic ray-depth and cell-cut oracles;
- `generic_ray_extension.m`, `generic_ray_depths.m`, `numerical_polar_cut.m`: black-box numerical cut path;
- `run_lifted_experiments.m`: numerical experiments and archived outputs.
- `run_jgo_solver_experiments.m`: independent local-solver and small-scale global checks.

## 7. Software requirements

- MATLAB with Optimization Toolbox (`linprog` and `fmincon`) and `convhulln` available;
- a MATLAB release supporting the language features used by the scripts;
- sufficient memory to materialize all rays only for the selected small or medium reference instances.

The large-dimensional generator and its analytic certificates do not require full ray materialization.
