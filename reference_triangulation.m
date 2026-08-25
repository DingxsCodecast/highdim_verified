function tri = reference_triangulation(n,q,varargin)
%REFERENCE_TRIANGULATION Validated small-scale Qhull triangulation of H_(n,q).

parser=inputParser;
addParameter(parser,'CoverageSamples',200,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&rem(x,1)==0);
addParameter(parser,'Seed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
parse(parser,varargin{:});
labels=enumerate_q_subsets(n,q);
H=zeros(size(labels,1),n-1);
for k=1:size(labels,1), H(k,:)=subset_to_h(labels(k,:),n).'; end

timer=tic;
cells=delaunayn(H,{'QJ','Qbb','Qc','Qs','Pp'});
qhull_seconds=toc(timer);
cells=sort(cells,2);
cells=unique(cells,'rows','stable');

rank_tolerances=zeros(size(cells,1),1);
sigma_min=zeros(size(cells,1),1);
volumes=zeros(size(cells,1),1);
valid=false(size(cells,1),1);
for k=1:size(cells,1)
    V=H(cells(k,:),:);
    B=(V(2:end,:)-V(1,:)).';
    s=svd(B);
    rank_tolerances(k)=max(size(B))*eps(norm(B,2));
    sigma_min(k)=s(end);
    valid(k)=sum(s>rank_tolerances(k))==n-1;
    if valid(k), volumes(k)=abs(det(B))/factorial(n-1); end
end
discarded=sum(~valid);
cells=cells(valid,:); sigma_min=sigma_min(valid); volumes=volumes(valid); rank_tolerances=rank_tolerances(valid);

analytic_volume=analytic_hypersimplex_volume(n,q);
volume_sum=sum(volumes);
volume_relative_error=abs(volume_sum-analytic_volume)/max(analytic_volume,realmin);

old_rng=rng; cleanup_rng=onCleanup(@()rng(old_rng));
rng(parser.Results.Seed,'twister');
P=zeros(parser.Results.CoverageSamples,n-1);
for k=1:size(P,1)
    count=min(10,size(H,1));
    idx=randperm(size(H,1),count);
    weights=-log(max(realmin,rand(count,1))); weights=weights/sum(weights);
    P(k,:)=weights.'*H(idx,:);
end
membership=tsearchn(H,cells,P);
coverage_rate=mean(~isnan(membership));

tri=struct();
tri.n=n; tri.q=q; tri.labels=labels; tri.H=H; tri.cells=cells;
tri.lifted_cell_labels=[cells,zeros(size(cells,1),1)]; % final zero denotes Z ray
tri.cell_volumes=volumes; tri.sigma_min=sigma_min; tri.rank_tolerances=rank_tolerances;
tri.number_cells=size(cells,1); tri.discarded_rank_deficient=discarded;
tri.analytic_volume=analytic_volume; tri.volume_sum=volume_sum;
tri.volume_relative_error=volume_relative_error; tri.coverage_rate=coverage_rate;
tri.uncovered_points=P(isnan(membership),:); tri.qhull_seconds=qhull_seconds;
tri.passed=discarded==0&&volume_relative_error<=1e-8&&coverage_rate==1;
if tri.passed
    tri.message='Reference triangulation passed rank, volume, and sampled coverage checks.';
else
    tri.message=sprintf('Triangulation failed: discarded=%d, volume relerr=%.3g, coverage=%.3g.', ...
        discarded,volume_relative_error,coverage_rate);
end
end
