setpaths;

myCluster  = parcluster();
parpool(myCluster);

for i=1:4
    switch i
        case 1
            n=500;p=50;
        case 2
            n=1000;p=100;
        case 3
            n=2000;p=200;
        case 4
            n=3000;p=300;               
    end
    d=2;
    clear beta
    beta(:,1)=[1,0,0,0,0,0,zeros(1,p-6)]';
    beta(:,2)=[0,1,0,0,0,0,zeros(1,p-6)]';
    time=zeros(20,2);
    angle=zeros(20,2);
    fun=zeros(20,2);
    parfor k=1:20
        x=randn(n,p);
        temp=x*beta;
        y=(temp(:,1)).^2+temp(:,2)+0.1*randn(n,1);
        initial=find_initial(x,y,d);
        
        tic;
        b0= dcsol2(x,y,d);
        t0=toc;
        
        [b2,~ ,t2]= MMRN(x,y,initial,struct('d',d,'verbosity',0,'tolnorm',1e-7,'maxiter',1000));
        
        fun(k,:)=[DistCov(x*b0,y),DistCov(x*b2,y)];
        time(k,:)=[t0,t2];
        angle(k,:)=[subspace(b0,beta),subspace(b2,beta)];        
    end
    save(['ModelA_Case_',num2str(i)],'time','fun','angle');    
end