function neighbors = baseline_nminus1_neighbors(S,n)
%BASELINE_NMINUS1_NEIGHBORS Deliberately incomplete simple-section diagnostic.
%
% This is not a Porembski/CDP implementation. It returns only the first n-1
% exact neighbors to quantify the consequence of hard-coding simplicity.
all_neighbors=delta_neighbors(S,n);
neighbors=all_neighbors(1:min(n-1,size(all_neighbors,1)),:);
end
