function labels = deterministic_subset_sample(n,q,requested)
%DETERMINISTIC_SUBSET_SAMPLE Reproducible q-subset labels without full enumeration.
count=min(requested,nchoosek(n,q));
if nchoosek(n,q)<=1e4
    all_labels=enumerate_q_subsets(n,q);
    idx=unique(round(linspace(1,size(all_labels,1),count)),'stable');
    labels=all_labels(idx,:);
else
    labels=zeros(count,q);
    labels(1,:)=1:q;
    for k=2:count
        labels(k,:)=sort(mod((0:q-1)+(k-1),n)+1);
    end
    labels=unique(labels,'rows','stable');
end
end
