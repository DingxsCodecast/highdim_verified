function value = relative_residual(actual, expected)
%RELATIVE_RESIDUAL Scale-aware Frobenius residual.
value = norm(actual-expected,'fro')/max(1,norm(expected,'fro'));
end
