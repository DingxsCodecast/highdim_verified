function value = lifted_reduced_value(x, instance, coordinates)
%LIFTED_REDUCED_VALUE Reduced value g in source or transformed coordinates.
if nargin<3, coordinates='source'; end
if strcmpi(coordinates,'transformed')
    x0=instance.generation.Mx\(x-instance.generation.shift_x);
    offset=instance.certificate.objective_offset;
else
    x0=x;
    offset=0;
end
n=instance.metadata.n;
ps=lifted_psi(x0(1:n-1),x0(n),instance.certificate);
value=x0(n)+x0(n+1)+instance.metadata.eta*ps+ ...
    min(0,3*(1-x0(n)-x0(n+1)))+offset;
end
