n=120;
p=24;
d=2;


Gamma(:,1)=[1,0,zeros(1,p-2)]';
Gamma(:,2)=[0,1,zeros(1,p-2)]';


beta=Gamma;


for i=1:p
    for j=1:p
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
    

x=randn(n,p)*spsig;
y=(x*beta(:,1))./(0.5+(x*beta(:,2)+1.5).^2)+0.2*randn(n,1);
%y=(x*beta(:,1))./(0.5+(x*beta(:,2)+1.5).^2)+0.2*trnd(8,n,1);


b0= dcsol2(x,y,d);
[b,st] = bic_sol(x,y,b0,0.5,200);

ang0(i)=subspace(beta,b0);
ang(i)=subspace(beta,b);

[fv1,beta1,st1] = cise(y,x,d,0.5,'SIR');
[fv2,beta2,st2] = cise(y,x,d,0.5,'PFC');

tpr=0;
tpr1=0;
tpr2=0;
for j=1:2
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
for j=3:p
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
 