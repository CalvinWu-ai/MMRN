n=60;
p=12;
d=2;


beta(:,1)=[0.5,0.5,0.5,0.5,zeros(1,p-4)]';  
beta(:,2)=[0.5,-0.5,0.5,-0.5,zeros(1,p-4)]';



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


parfor i=1:120
    

x=randn(n,p)*spsig;

y1=(x*beta(:,1)) +randn(n,1);
y2=(x*beta(:,2)+0.5).^2 + randn(n,1);

%y1=(x*beta(:,1)) +trnd(8,n,1);
%y2=(x*beta(:,2)+0.5).^2 + trnd(8,n,1);

y=[y1,y2];


b0 = PRSIR1(y,x,d,6);

[b1,st1] = bic_sol(x,y,b0,0.5,100);
[b2,st2] = aic_sol(x,y,b0,0.5,100);

tpr1=0;
tpr2=0;

for j=1:4
      if (st1(j)==1)
         tpr1=tpr1+1;
      end 
      
       if (st2(j)==1)
         tpr2=tpr2+1;
      end 
      
end

fpr1=0;
fpr2=0;

for j=5:p
      if (st1(j)==0)
         fpr1=fpr1+1;
      end
        if (st2(j)==0)
         fpr2=fpr2+1;
      end
      
      
end
 
 stpr1(i)=tpr1;
 sfpr1(i)=fpr1;
 
 stpr2(i)=tpr2;
 sfpr2(i)=fpr2;
 
end

delete(gcp);
 