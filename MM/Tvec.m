function T=Tvec(n)

% Compute T_n s.t. vec(W')=T_n vec(W), W \in R^{n*n}

coord=zeros(n^2,3);
for k=1:n^2
    r=mod(k,n);
    p=(k-r)/n;
    if r==0
        row=p;
        col=n;
    else
        row=p+1;
        col=r;
    end
    
    index=(col-1)*n+row;
    coord(k,:)=[k,index,1];      
end
T = sparse(coord(:,1),coord(:,2),coord(:,3),n^2,n^2);
end