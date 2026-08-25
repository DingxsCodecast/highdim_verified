function h = subset_to_h(S, n)
%SUBSET_TO_H Project a hypersimplex subset label by deleting coordinate n.

S = sort(S(:).');
if any(S < 1) || any(S > n) || numel(unique(S)) ~= numel(S)
    error('subset_to_h:InvalidSubset', 'S must contain distinct indices in 1:n.');
end
h = zeros(n-1,1);
h(S(S < n)) = 1;
end
