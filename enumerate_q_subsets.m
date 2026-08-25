function subsets = enumerate_q_subsets(n, q)
%ENUMERATE_Q_SUBSETS Lexicographically enumerate all q-subsets of 1:n.

validateattributes(n, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(q, {'numeric'}, {'scalar','integer','>=',0,'<=',n});
subsets = nchoosek(1:n, q);
end
