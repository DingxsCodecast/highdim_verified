function result = hrep_adjacency_graph(A, rhs, vertices, varargin)
%HREP_ADJACENCY_GRAPH Recover vertex adjacency using common active-normal rank.
%   Two distinct vertices are adjacent exactly when the common active
%   inequalities define a one-dimensional minimal face.

parser = inputParser;
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(parser,varargin{:});
tol = parser.Results.Tolerance;

validateattributes(A,{'numeric'},{'2d','real','finite'});
validateattributes(rhs,{'numeric'},{'column','real','finite','numel',size(A,1)});
validateattributes(vertices,{'numeric'},{'2d','real','finite'});
if size(vertices,1) ~= size(A,2)
    error('hrep_adjacency_graph:DimensionMismatch', ...
        'Each vertex column must have size(A,2) entries.');
end
d = size(A,2);
nv = size(vertices,2);
residual = A*vertices-rhs;
scale = 1 + abs(rhs) + abs(A)*abs(vertices);
active = abs(residual) <= tol*scale;
adjacency = false(nv,nv);
face_dimension = nan(nv,nv);
common_active_rank = nan(nv,nv);
for i = 1:nv
    face_dimension(i,i) = 0;
    common_active_rank(i,i) = d;
    for j = i+1:nv
        rows = active(:,i) & active(:,j);
        if any(rows)
            common = A(rows,:);
            s = svd(common);
            rank_tol = max(size(common))*eps(max(1,norm(common,2)));
            common_rank = sum(s > rank_tol);
        else
            common_rank = 0;
        end
        dimension = d-common_rank;
        common_active_rank(i,j) = common_rank;
        common_active_rank(j,i) = common_rank;
        face_dimension(i,j) = dimension;
        face_dimension(j,i) = dimension;
        adjacency(i,j) = dimension == 1;
        adjacency(j,i) = adjacency(i,j);
    end
end

result = struct('adjacency',adjacency,'degree',sum(adjacency,2), ...
    'active',active,'active_count',sum(active,1).', ...
    'face_dimension',face_dimension,'common_active_rank',common_active_rank, ...
    'dimension',d,'number_vertices',nv,'tolerance',tol);
end
