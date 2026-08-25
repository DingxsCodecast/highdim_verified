function minima = lifted_row_minima(M, n, q)
%LIFTED_ROW_MINIMA Row minima over X_(0,q) times [0,1].

if size(M,2) ~= n+1
    error('lifted_row_minima:DimensionMismatch', 'M must have n+1 columns.');
end
minima = pyramid_row_minima(M(:,1:n), n, q) + min(0, M(:,n+1));
end
