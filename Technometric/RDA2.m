RawX   = load('rawX.csv');
Y=[143,84,98,83,153,141,191,130,744,381,1047,806,621,849,475,966,708,487,1447,1693,1731,1025,376,126,102,149,91,153,235,68];
RawX=RawX';
Y=Y';

mX = bsxfun(@minus, RawX, mean(RawX));
sX= bsxfun(@rdivide, mX, std(mX));

H=6;
m=30;
rang=2;
setpaths

[nx1,index1] = dc_ordering(Y,sX);
[nx2,index2]=seq_dc1(nx1,Y,m,H,rang);
idx=index1(index2);

b0= dcsol2(nx2,Y,1);
b0=-b0;
[b,st] = aic_sol(nx2,Y,b0,rang,500);

L=logical(st);
nx3=nx2(:,L);
sdx=idx(L); % the selected index

beta=dcsol2(nx3,Y,1); % the estimated coefficients of selected variables
