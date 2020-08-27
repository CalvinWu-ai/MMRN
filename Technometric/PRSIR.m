function [b0,d,G,N] = PRSIR(y,x,H)

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

chicalc=zeros(H,1);
chitab=zeros(H,1);

for i= 0:(H-2)
	chicalc(i+1)=n*sum(E((i+1):p));
	df=(H-i-1)*(p-i);
	chitab(i+1)=chi2inv(0.95,df);
end

	d=0;
	k=1;
    
    while (chicalc(k)>chitab(k) && d < H-1)
       d=d+1;
       k=k+1;
    end
    
    if (d==0)
        b0=1;
    else
	b0= iqN*V(:,1:d);
    end
    

