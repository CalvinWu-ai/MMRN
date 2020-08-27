function b= dcsol1(x,y,b0)

% compute the original solution without penalty function
[n,p]=size(x);
X=x-ones(n,1)*mean(x);
N=cov(X);
N1=sqrtm(N);
N2=inv(N1);

v0=N1*b0;

options = optimset('Display','notify','Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-4,'TolX',1e-4,'LargeScale', 'off');

v = fmincon(@(c)dcreg(c,x,y,N2),v0,[],[],[],[],[],[],@tcon,options);

b= N2*v;

return;