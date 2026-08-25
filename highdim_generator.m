function [Q, A, B, rhs_A, rhs_B, c, d, minimum, solution_x, solution_y, certificate, instance] = highdim_generator(n, varargin)
%HIGHDIM_GENERATOR Positional-output wrapper for existing DBLP solver code.
%
%   The first ten outputs follow smart_generator's convention.  The last
%   two outputs expose the certificate and the complete instance struct.

instance = generate_highdim_dblp(n, varargin{:});
Q = instance.Q;
A = instance.A;
B = instance.B;
rhs_A = instance.rhs_A;
rhs_B = instance.rhs_B;
c = instance.c;
d = instance.d;
minimum = instance.minimum;
solution_x = instance.solution_x;
solution_y = instance.solution_y;
certificate = instance.certificate;
end
