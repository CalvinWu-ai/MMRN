function [nx,index] = dc_ordering(y,x)

[n,p]=size(x);

for i=1:p
    [temp,dc(i)]=DistCorrVec(x(:,i),y);
end

[a,index]=sort(dc);

nx=x(:,index);
