function b0 = PRSIR1(y,x,d,H)

m=1000;

[n,q]=size(y);
[n,p]=size(x);

X=x-ones(n,1)*mean(x);
N=cov(X);
iqN=inv(sqrtm(N));
sX= X * iqN;


G=zeros(p,p);

for i=1:m
    
    t=orth(randn(q,1));
    
    G=G+SIRm(y*t,sX,H);
end

G=G/m;

[V,D]=eig(G);
[E,I] = sort(diag(D),'descend');
V = V(:, I);

b0= iqN*V(:,1:d);