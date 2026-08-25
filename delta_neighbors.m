function neighbors = delta_neighbors(S, n)
%DELTA_NEIGHBORS Exact exchange-one-index neighbors in Delta(n,q).

S = sort(S(:).');
q = numel(S);
if q < 1 || q >= n || any(S < 1) || any(S > n) || numel(unique(S)) ~= q
    error('delta_neighbors:InvalidSubset', 'S must be a proper subset of 1:n.');
end
outside = setdiff(1:n, S, 'stable');
neighbors = zeros(q*(n-q), q);
cursor = 0;
for i = 1:q
    for j = 1:numel(outside)
        cursor = cursor + 1;
        T = S;
        T(i) = outside(j);
        neighbors(cursor,:) = sort(T);
    end
end
neighbors = sortrows(neighbors);
end
