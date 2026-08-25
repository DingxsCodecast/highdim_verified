function value = lifted_psi(u,t,certificate)
%LIFTED_PSI Evaluate the normalized exposing functional.
value=(certificate.w(:)'*u(:)-certificate.m*t)/certificate.D;
end
