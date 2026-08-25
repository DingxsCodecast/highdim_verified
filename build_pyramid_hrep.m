function [A, rhs] = build_pyramid_hrep(n, q)
%BUILD_PYRAMID_HREP H-representation of the pyramid over Delta(n,q).

validateattributes(n, {'numeric'}, {'scalar','integer','>=',4});
validateattributes(q, {'numeric'}, {'scalar','integer','>=',2,'<=',n-2});
d = n - 1;
A = [
    -eye(d), zeros(d,1);
     eye(d), -ones(d,1);
     ones(1,d), -q;
    -ones(1,d), q-1;
     zeros(1,d), 1
];
rhs = [zeros(2*n,1); 1];
end
