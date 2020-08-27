n=60;
p=24;
d=2;

Gamma(:,1)=[0.5,0.5,0.5,0.5,zeros(1,p-4)]';
Gamma(:,2)=[0.5,-0.5,0.5,-0.5,zeros(1,p-4)]';

beta=Gamma;



for i=1:p-1
    for j=1:p-1
        psig(i,j)=0.5^abs(i-j);
    end
end

setpaths;
spsig = sqrtm(psig);
myCluster = parcluster('local');
myCluster.NumWorkers = 30;
setpaths;
parpool('local',30);


parfor i=1:400
 
x210 =randn(n,p-1)*spsig;
x01=abs(x210(:,1)+x210(:,2))+randn(n,1);
x=[x01,x210];

y=(x*beta(:,1)).^2+abs(x*beta(:,2))+0.5*randn(n,1);
%y=(x*beta(:,1)).^2+abs(x*beta(:,2))+0.5*trnd(8,n,1);

b0= dcsol2(x,y,d);
[b,st] = bic_sol(x,y,b0,0.5,200);

[fv1,beta1,st1] = cise(y,x,d,0.5,'SIR');

[fv2,beta2,st2] = cise(y,x,d,0.5,'PFC');


tpr=0;
tpr1=0;
tpr2=0;
for j=1:4
      if (st(j)==1)
         tpr=tpr+1;
      end
      if (st1(j)==1)
         tpr1=tpr1+1;
      end      
      if (st2(j)==1)
         tpr2=tpr2+1;
      end       
      
end

fpr=0;
fpr1=0;
fpr2=0;
for j=5:p
      if (st(j)==0)
         fpr=fpr+1;
      end
      if (st1(j)==0)
         fpr1=fpr1+1;
      end
      if (st2(j)==0)
         fpr2=fpr2+1;
      end      
      
      
end
 
 stpr(i)=tpr;
 sfpr(i)=fpr;
 
  stpr1(i)=tpr1;
 sfpr1(i)=fpr1;
 
  stpr2(i)=tpr2;
 sfpr2(i)=fpr2;
 
end

delete(gcp);
 