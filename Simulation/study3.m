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

myCluster  = parcluster();
parpool(myCluster);



for n=[60 120]
    time=zeros(100,2);
    tpr=zeros(100,2); 
    fpr=zeros(100,2);
    angle=zeros(100,3);
    parfor i=1:100        
        x210 =randn(n,p-1)*spsig;
        x01=abs(x210(:,1)+x210(:,2))+randn(n,1);
        x=[x01,x210];
        temp=x*beta;
        y=(temp(:,1)).^2+abs(temp(:,2))+0.5*randn(n,1);
        
        initial=find_initial(x,y,d);
        b0= MMRN(x,y,initial,struct('d',d,'verbosity',0,'tolnorm',1e-7,'maxiter',1000));
        
        tic;
        [b1,st1] = bic_sol(x,y,b0,0.5,200);
        t1=toc;
        
        [b2,st2,t2] = PMMRNBic(x,y,b0,0.5);
        
      
    
        time(i,:)=[t1,t2];
        angle(i,:)=[subspace(b0,beta),subspace(b1,beta),subspace(b2,beta)];

        tpr1=sum(st1(1:4))/4;
        tpr2=sum(st2(1:4))/4;
        tpr(i,:)=[tpr1,tpr2];

        fpr1=sum(st1(5:p))/(p-4);
        fpr2=sum(st2(5:p))/(p-4);
        fpr(i,:)=[fpr1,fpr2];
    end
    save(['Study3_',num2str(n)],'time','tpr','fpr','angle');      
end

delete(gcp);
