function D=veck_to_vec(n)

% Compute D_n s.t. vec(U)=D_n veck(U), U \in Skew(n)

coord=zeros(n^2,3);
for k=1:n^2
    r=mod(k,n);
    p=(k-r)/n;
    if r==0
        row=n;
        col=p;
    else
        row=r;
        col=p+1;
    end
    
    if row>col
        index=(col-1)*n+row-col*(col+1)/2;
        num=1;
    elseif row<col
        temp=row;
        row=col;
        col=temp;        
        index=(col-1)*n+row-col*(col+1)/2;
        num=-1;
    else
        index=0;
        num=0;
    end
    coord(k,:)=[k,index,num];    
end

coord(coord(:,2)==0,:)=[];
D = sparse(coord(:,1),coord(:,2),coord(:,3),n^2,n*(n-1)/2);

end

