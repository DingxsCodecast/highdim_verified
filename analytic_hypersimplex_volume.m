function volume = analytic_hypersimplex_volume(n,q)
%ANALYTIC_HYPERSIMPLEX_VOLUME Euclidean volume in projected coordinates.
volume=eulerian_number(n-1,q-1)/factorial(n-1);
end
