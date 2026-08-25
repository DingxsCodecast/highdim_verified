function [A, rhs] = build_lifted_hrep(n, q)
%BUILD_LIFTED_HREP H-representation of X_(0,q) times [0,1].

[A0, rhs0] = build_pyramid_hrep(n, q);
A = [
    A0, zeros(size(A0,1),1);
    zeros(1,n), -1;
    zeros(1,n), 1
];
rhs = [rhs0; 0; 1];
end
