function M = make_affine_map(dimension, orientation, reflections, scaling_range)
%MAKE_AFFINE_MAP Orthogonal mixing followed by positive diagonal scaling.

validateattributes(dimension, {'numeric'}, {'scalar','integer','>=',1});
orientation = upper(char(orientation));
switch orientation
    case 'NONE'
        U = eye(dimension);
        M = U;
        return;
    case 'HH'
        validateattributes(reflections, {'numeric'}, {'scalar','integer','>=',1});
        U = eye(dimension);
        for k = 1:reflections
            v = randn(dimension,1);
            v = v / norm(v,2);
            U = (eye(dimension) - 2*(v*v')) * U;
        end
    case 'QR'
        [U,R] = qr(randn(dimension),0);
        signs = sign(diag(R));
        signs(signs == 0) = 1;
        U = U * diag(signs);
    otherwise
        error('make_affine_map:InvalidOrientation', ...
            'Orientation must be ''none'', ''HH'', or ''QR''.');
end
scales = scaling_range(1) + diff(scaling_range)*rand(dimension,1);
M = U * diag(scales);
end
