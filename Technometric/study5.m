n=200;
p=1000;  


for i=1:p
    for j=1:p
        psig(i,j)=0.5^abs(i-j);
    end
end


H=6;
m=30;

%setting

setpaths;
parpool('local',30);
spsig = sqrtm(psig);
spsig=(spsig+spsig')/2;
myCluster = parcluster('local');
myCluster.NumWorkers = 30;
setpaths;
parpool('local',30);


parfor i=1:100

rng('shuffle');
x=randn(n,p)*spsig;

y = 1+ exp(x(:,501)+x(:,502)+x(:,503)+x(:,504)) + randn(n,1);

[nx1,index1] = dc_ordering(y,x);
[nx2,index2]=seq_dc1(nx1,y,m,H,2);
idx=index1(index2);

in=length(idx);

it=0;
for j=1:in
    if(idx(j)==1 || idx(j)==2 || idx(j)==3 || idx(j)==4 || idx(j)==7 ||idx(j)==8||idx(j)==9||idx(j)==10)
         it=it+1;
    end
end
    
tt(i)=it;

b0= dcsol2(nx2,y,2);
[beta,st] = aic_sol(nx2,y,b0,0.3,100);


fs=length(st);
tpr=0;
for j=1:fs
      if (st(j)==1 && (idx(j)==1 || idx(j)==2 || idx(j)==3 || idx(j)==4 || idx(j)==7 ||idx(j)==8||idx(j)==9||idx(j)==10))
         tpr=tpr+1;
      end    
end

tr(i)=in;
stpr(i)=tpr;
sfpr(i)=sum(st)-tpr;
 
end

delete(gcp);

