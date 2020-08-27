n=200;
p=1000;  
d=2;

for i=1:p
    for j=1:p
        psig(i,j)=0.5^abs(i-j);
    end
end


H=6;
m=30;

%setting

setpaths;
myCluster = parcluster('local');
myCluster.NumWorkers = 30;
parpool('local',30);

spsig = sqrtm(psig);
spsig=(spsig+spsig')/2;



parfor i=1:100

rng('shuffle');
x=randn(n,p)*spsig;


y=0.4*(0.5*x(:,1)+0.5*x(:,2)+0.5*x(:,3)+0.5*x(:,4)).^2 + 3*sin(0.25*x(:,1)-0.25*x(:,2)+0.25*x(:,3)-0.25*x(:,4))+0.2*randn(n,1);


[nx1,index1] = dc_ordering(y,x);
[nx2,index2]=seq_dc1(nx1,y,m,H);
idx=index1(index2);

in=length(idx);

it=0;
for j=1:in
    if (idx(j)==1 || idx(j)==2 || idx(j)==3 || idx(j)==4)
         it=it+1;
      end
end
    
tt(i)=it;

b0= dcsol2(nx2,y,d);
[beta,st] = aic_sol(nx2,y,b0,0.3,100);


fs=length(st);
tpr=0;
for j=1:fs
      if (st(j)==1 && (idx(j)==1 || idx(j)==2 || idx(j)==3 || idx(j)==4))
         tpr=tpr+1;
      end    
end

tr(i)=in;
stpr(i)=tpr;
sfpr(i)=sum(st)-tpr;
 
end

delete(gcp);

