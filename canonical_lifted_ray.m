function [source_ray, transformed_ray] = canonical_lifted_ray(instance, ray_label)
%CANONICAL_LIFTED_RAY Return unnormalized adjacent-vertex difference.
n=instance.metadata.n;
if ischar(ray_label)||isstring(ray_label)
    if ~strcmpi(char(ray_label),'Z')
        error('canonical_lifted_ray:InvalidLabel','String ray label must be Z.');
    end
    source_ray=[zeros(n,1);1];
else
    source_ray=[subset_to_h(ray_label,n);1;0];
end
transformed_ray=instance.generation.Mx*source_ray;
end
