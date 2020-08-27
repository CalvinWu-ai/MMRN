n=120;
p=24;
d=1;


beta=[1,1,1,1,zeros(1,20)]'/sqrt(4); 


for i=1:p
    for j=1:p
        psig(i,j)=0.5^abs(i-j);
    end
end

spsig = sqrtm(psig);

setpaths;
parpool('local',30);


parfor i=1:400
    

x=randn(n,p)*spsig;


y=(x*beta+0.5).^2+0.5*randn(n,1);

b0= dcsol2(x,y,2);
[b2,st] = bic_sol(x,y,b0,0.5,200);

[fv1,beta1,st1] = cise(y,x,1,0.5,'SIR');
[fv2,beta2,st2] = cise(y,x,1,0.5,'PFC');


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
 