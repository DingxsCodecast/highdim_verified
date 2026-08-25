# 可验证的高维退化 DBLP 基准（structural + lifted）

本目录与论文 `highdim_dblp_generator_lifted_complete.tex` 的符号和坐标约定一致。它提供两个共享 `Delta(n,q)` 几何的模式。

- Structural：`X_(0,q)` 是 `Delta(n,q)` 的棱锥，退化 apex 是唯一全局最优解，适合邻接与几何验证。
- Lifted：`Xtilde_(0,q)=X_(0,q) x [0,1]`，lifted apex 是严格、非全局 PGM；另有已知唯一全局最优解、指定值 `alpha` 的可行 incumbent 和每条 incident ray 的闭式最大 multiplier，适合 decomposition / polar-cut 实验。

合法参数为 `n>=4`、`2<=q<=n-2`。lifted 模式还要求 `0<eta<2`、`-1<alpha<0`。

## 1. Structural 模式

```matlab
addpath('E:\Seafile\Warehouse\Submission\WaitingList\Generator\cpm_codes\generator\highdim_verified');

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

兼容原求解器位置输出约定：

```matlab
[Q,A,B,rhs_A,rhs_B,c,d,fstar,xstar,ystar,certificate,instance] = ...
    highdim_generator(20,'Q',3,'Seed',20260820);
```

## 2. Lifted 模式

```matlab
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

在 `Orientation='none'` 时不应用缩放、旋转或平移，因此 source 值严格对应论文的 `0, alpha, -1`。HH/QR 模式中目标函数整体增加保存于 `L.certificate.objective_offset` 的常数；三者的差值和 canonical ray multiplier 不变。

canonical ray 使用“PGM 到相邻顶点的差向量”，不做单位长度归一化：

```matlab
Slabel = L.certificate.Sstar;
[r0,rT] = canonical_lifted_ray(L,Slabel);
lambdaS = lambda_oracle(L,Slabel);
lambdaZ = lambda_oracle(L,'Z');
```

## 3. 独立邻接、三角剖分和 cell cut

```matlab
N = delta_neighbors([1 3 5],8);       % exact exchange-one-index oracle
B = baseline_nminus1_neighbors([1 3 5],8); % 故意不完整的 n-1 diagnostic；不是 CDP

% 不使用 subset/exchange 规则：从 H-description 枚举顶点，
% 再用共同 active normals 的秩判定最小公共 face 的维数。
[AH,bH] = hypersimplex_hrep(8,3);
[VH,enumInfo] = enumerate_hrep_vertices(AH,bH);
GH = hrep_adjacency_graph(AH,bH,VH);
assert(enumInfo.number_vertices == nchoosek(8,3))
assert(all(GH.degree == 3*(8-3)))

% Delta(4,2) 上的过程级 neighbor-count witness。
W = porembski_first_step_witness();
assert(W.passed && W.actual_neighbors == 4 && W.n_minus_one_neighbors == 3)

T = reference_triangulation(8,3, ...
    'CoverageSamples',500,'Seed',20260820);
assert(T.passed,T.message)

L0 = generate_lifted_highdim_dblp(8,'Q',3, ...
    'Eta',0.5,'Alpha',-0.5,'Orientation','none');
C = oracle_polar_cut(L0,T,1,'Coordinates','source');
assert(C.passed)
disp(C.lambda_oracle)
disp(C.lambda_from_cut)

% 独立主实验路径：将 reduced value 当作黑箱，用 bracket+bisection 求深度；
% numerical_polar_cut 不调用 lambda_oracle。
D = generic_ray_depths(L0,T.labels,'Coordinates','source');
CN = numerical_polar_cut(L0,T,1,D,'Coordinates','source');
assert(D.passed && CN.passed)
```

`reference_triangulation` 仅用于小规模：它物化全部 `nchoosek(n,q)` 条 hypersimplex rays，用 Qhull `QJ` 构造候选剖分，再独立检查每个 cell 的秩、Eulerian-number 解析体积和随机内部点覆盖。大规模实例应使用 ray sampling 或隐式算法。

cell cut 是 branch-specific 的：`a'*(x-xbar)>=1` 必须和当前 simplicial cell 约束一起使用，不能把一个 cell 的 cut 不加区分地施加到其它 cell。

## 4. 自动验证与测试

```powershell
matlab -batch "addpath('E:\Seafile\Warehouse\Submission\WaitingList\Generator\cpm_codes\generator\highdim_verified'); run_generator_tests"
matlab -batch "addpath('E:\Seafile\Warehouse\Submission\WaitingList\Generator\cpm_codes\generator\highdim_verified'); stress_test_highdim_generator"
```

回归测试包含 15 项：一般 `q`、HH/QR、手工 smoke identities、exchange
oracle、独立 H-description 顶点/邻接恢复、四维过程级 witness、黑箱射线延拓、
解析体积覆盖、source/transformed numerical cut、RNG 可重复性、旧包装接口和
非法输入。压力测试交叉 12 组 `(n,q)`、3 个 seed、2 种 orientation 和 2 个
mode，共 144 个实例。

## 5. 复现实验

```matlab
out = fullfile(pwd,'results');
run_lifted_experiments(out,'Profile','full');
```

输出：

- `experiment_results.mat`：全部 MATLAB 结构和 raw metrics；
- `table_neighborhood.csv`：exact oracle 与 `n-1` diagnostic；
- `table_independent_neighborhood.csv`：从 H-description 独立恢复的顶点和边；
- `table_procedure_witness.csv`：`Delta(4,2)` 上的实际 neighbor-count witness；
- `table_decomposition.csv`：cell 数、秩、解析体积和覆盖；
- `table_cut_depth.csv`：96 个组合的黑箱 bisection depth、解析误差和 numerical cut 有效性；
- `table_orientation.csv`：none/HH/QR、3 seeds 的 affine robustness；
- `table_generation.csv`：两类实例的 5-seed 生成加验证时间；
- `table_pipeline.csv`：确定性 search-to-cut 流程；
- `environment.txt`：CPU、内存、Windows、MATLAB、并行池和 GPU 状态。

正式实验容差为 `tol_feas=tol_cut=tol_cone=1e-9`；rank tolerance 使用 `max(size(R))*eps(norm(R,2))`。任何失败都应保留 subset/cell labels、cut 系数和相关残差，不能只记录 `status=fail`。

## 6. 主要文件

- `generate_highdim_dblp.m` / `validate_highdim_dblp.m`：structural 模式；
- `generate_lifted_highdim_dblp.m` / `validate_lifted_highdim_dblp.m`：lifted 模式；
- `delta_neighbors.m`：邻接 oracle；
- `hypersimplex_hrep.m` / `enumerate_hrep_vertices.m` /
  `hrep_adjacency_graph.m`：独立 H-description 邻接验证；
- `porembski_first_step_witness.m`：四维过程级 neighbor-count witness；
- `reference_triangulation.m`：小规模 reference decomposition；
- `lambda_oracle.m` / `oracle_polar_cut.m`：ray-depth 与 cell cut；
- `generic_ray_extension.m` / `generic_ray_depths.m` /
  `numerical_polar_cut.m`：不使用闭式深度的黑箱主实验路径；
- `run_lifted_experiments.m`：论文数值实验与归档。
