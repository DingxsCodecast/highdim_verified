function value = eulerian_number(m,k)
%EULERIAN_NUMBER Eulerian number A(m,k) using the standard recurrence.
validateattributes(m,{'numeric'},{'scalar','integer','>=',1});
validateattributes(k,{'numeric'},{'scalar','integer','>=',0});
if k>=m, value=0; return; end
A=zeros(m,m);
A(1,1)=1;
for i=2:m
    for j=0:i-1
        left=0; right=0;
        if j>0, left=(i-j)*A(i-1,j); end
        if j<=i-2, right=(j+1)*A(i-1,j+1); end
        A(i,j+1)=left+right;
    end
end
value=A(m,k+1);
end
