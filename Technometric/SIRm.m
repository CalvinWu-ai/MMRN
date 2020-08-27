function M = SIRm(y,x,nslice)

[n,p]=size(x);

[a, index] = sort(y);
xmean = mean(x);
group = sort(mod(1:n, nslice)); 

sigmaeta = zeros(p, p);
diffmean = zeros(1, p);

for i=0:(nslice-1)
    diffmean = mean( x(index(group==i), :) ) - xmean;
    sigmaeta = sigmaeta + mean(group == i) * diffmean' * diffmean;
end

M = sigmaeta;
