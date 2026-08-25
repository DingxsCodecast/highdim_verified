function [vertices, info] = enumerate_hrep_vertices(A, rhs, varargin)
%ENUMERATE_HREP_VERTICES Generic small-scale vertex enumeration from Ax <= rhs.
%   The routine enumerates full-rank d-row bases, solves their intersections,
%   retains feasible points, and merges duplicates within a numerical tolerance.

parser = inputParser;
addParameter(parser,'Tolerance',1e-9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(parser,varargin{:});
tol = parser.Results.Tolerance;

validateattributes(A,{'numeric'},{'2d','real','finite'});
validateattributes(rhs,{'numeric'},{'column','real','finite','numel',size(A,1)});
[m,d] = size(A);
if m < d
    error('enumerate_hrep_vertices:TooFewRows', ...
        'At least d inequality rows are required in dimension d.');
end

bases = nchoosek(1:m,d);
candidates = zeros(d,size(bases,1));
accepted = 0;
full_rank_bases = 0;
for k = 1:size(bases,1)
    AI = A(bases(k,:),:);
    s = svd(AI);
    rank_tol = max(size(AI))*eps(max(1,norm(AI,2)));
    if sum(s > rank_tol) ~= d
        continue;
    end
    full_rank_bases = full_rank_bases + 1;
    x = AI\rhs(bases(k,:));
    feasibility_scale = 1 + abs(rhs) + sum(abs(A),2)*max(1,norm(x,inf));
    if all(A*x-rhs <= tol*feasibility_scale)
        accepted = accepted + 1;
        candidates(:,accepted) = x;
    end
end
candidates = candidates(:,1:accepted);
if isempty(candidates)
    vertices = zeros(d,0);
else
    rows = uniquetol(candidates.',10*tol,'ByRows',true,'DataScale',1);
    rows(abs(rows) <= 10*tol) = 0;
    rows(abs(rows-1) <= 10*tol) = 1;
    rows = sortrows(rows,1:d);
    vertices = rows.';
end

max_violation = 0;
if ~isempty(vertices)
    max_violation = max(max(A*vertices-rhs,[],1));
end
info = struct('dimension',d,'inequalities',m,'candidate_bases',size(bases,1), ...
    'full_rank_bases',full_rank_bases,'accepted_bases',accepted, ...
    'number_vertices',size(vertices,2),'maximum_violation',max_violation, ...
    'tolerance',tol);
end
