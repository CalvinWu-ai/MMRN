function b= dcsol(x,y,d)

% compute the original solution without penalty function
[n,p]=size(x);
X=x-ones(n,1)*mean(x);
N=cov(X);
N1=sqrtm(N);
N2=inv(N1);
sX= X * N2;


options = optimset('Display','notify','Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-4,'TolX',1e-4,'LargeScale', 'off');

[WX,W1,auxmtx] = SIR(y,sX,'cont',d,'nslices',6);
v1=orth(W1);
[WX,W2,auxmtx] = DR(y,sX,'cont',d,'nslices',6);
v2=orth(W2);
W3 = dMAVE(sX, y, d);
v3=orth(W3);


p1=-DistCorrVec(x*N2*v1,y);
p2=-DistCorrVec(x*N2*v2,y);
p3=-DistCorrVec(x*N2*v3,y);


[C,I]=min([p1,p2,p3]);

if (I(1)==1)
    vi=v1;
elseif (I(1)==2)
    vi=v2;
else
    vi=v3;
end


S=zeros(p,d,500);

vc=orthcomp(vi);

for k=1:500
    
    t=orth(vi+0.1*vc*randn(p-d,d));
    S(:,:,k)=t;
    ts(k)=dcreg(t,x,y,N2);
end

[D,J]=min(ts);
v0=S(:,:,J(1));


v = fmincon(@(c)dcreg(c,x,y,N2),v0,[],[],[],[],[],[],@tcon,options);

b= N2*v;

return;