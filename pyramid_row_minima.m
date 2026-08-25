function minima = pyramid_row_minima(M, n, q)
%PYRAMID_ROW_MINIMA Row minima of M*x over the Delta(n,q) pyramid.

if size(M,2) ~= n
    error('pyramid_row_minima:DimensionMismatch', 'M must have n columns.');
end
minima = zeros(size(M,1),1);
for row = 1:size(M,1)
    augmented = [M(row,1:n-1), 0];
    ordered = sort(augmented, 'ascend');
    base_value = M(row,n) + sum(ordered(1:q));
    minima(row) = min(0, base_value);
end
end
