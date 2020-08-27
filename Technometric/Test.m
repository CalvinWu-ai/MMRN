r=6;
n=200;
p=1;
u=2;
eta=[2,1]';

setpaths;

Gamma(:,1)=ones(r,1);
c=[1,-1];
Gamma(:,2) = kron(ones(1,r/2),c)';

Gamma=Gamma/sqrt(r);

Gamma0=orth(orthcomp(Gamma));


om=1;
om0=3;

for k=1:400
    


x=randn(n,p);

y = x*eta'*Gamma' + randn(n,u)*om*Gamma'+randn(n,r-u)*om0*Gamma0';


[GX,b1,L,dhat]=mlm_fit(y,x,u);

dis1(k)=norm(b1*b1'-Gamma*Gamma','fro');



end