function [V, labels] = materialize_pyramid_vertices(n, q)
%MATERIALIZE_PYRAMID_VERTICES Columns are apex followed by labeled base vertices.

labels = enumerate_q_subsets(n,q);
V = zeros(n,size(labels,1)+1);
for k = 1:size(labels,1)
    V(:,k+1) = [subset_to_h(labels(k,:),n);1];
end
end
