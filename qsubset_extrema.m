function [minimum_value, maximum_value] = qsubset_extrema(weights, q)
%QSUBSET_EXTREMA Minimum and maximum q-subset sums, column by column.
%
% WEIGHTS is n-by-p. Outputs are 1-by-p.

validateattributes(weights, {'numeric'}, {'2d','finite'});
n = size(weights,1);
validateattributes(q, {'numeric'}, {'scalar','integer','>=',1,'<=',n});
ordered = sort(weights, 1, 'ascend');
minimum_value = sum(ordered(1:q,:), 1);
maximum_value = sum(ordered(end-q+1:end,:), 1);
end
