function [A, rhs] = hypersimplex_hrep(n, q)
%HYPERSIMPLEX_HREP Full-dimensional H-representation of projected Delta(n,q).

validateattributes(n, {'numeric'}, {'scalar','integer','>=',4});
validateattributes(q, {'numeric'}, {'scalar','integer','>=',2,'<=',n-2});
d = n - 1;
A = [
    -eye(d);
     eye(d);
     ones(1,d);
    -ones(1,d)
];
rhs = [zeros(d,1); ones(d,1); q; -(q-1)];
end
