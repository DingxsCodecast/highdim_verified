function lambda = lambda_oracle(instance, ray_label)
%LAMBDA_ORACLE Exact canonical multiplier for a subset ray or 'Z'.
alpha=instance.metadata.alpha_source;
if ischar(ray_label)||isstring(ray_label)
    if ~strcmpi(char(ray_label),'Z')
        error('lambda_oracle:InvalidLabel','String ray label must be Z.');
    end
    lambda=(3-alpha)/2;
else
    h=subset_to_h(ray_label,instance.metadata.n);
    psiS=lifted_psi(h,1,instance.certificate);
    lambda=(3-alpha)/(2-instance.metadata.eta*psiS);
end
end
