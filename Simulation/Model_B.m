setpaths;

myCluster  = parcluster();
parpool(myCluster);


for i=1:2
    for j=1:3
        if i==1
            n=100;
            p=6;
        else
            n=500;
            p=20;
        end
        
        clear beta
        beta(:,1)=[1,0,0,0,0,0,zeros(1,p-6)]';
        beta(:,2)=[0,1,0,0,0,0,zeros(1,p-6)]';
        d=2;
                
        time=zeros(100,2);
        angle=zeros(100,2);
        fun=zeros(100,2);
        parfor k=1:100
            if j==1
                x=randn(n,p);
            elseif j==2
                x=unifrnd(-2,2,[n,p])
            else
                x=binornd(10,0.1,[n,p]);           
            end
            temp=x*beta;
            y=sign(2*temp(:,1)+randn(n,1)).*log(abs(2*temp(:,2)+4+randn(n,1)));
            initial=find_initial(x,y,d);
            
            tic;
            b0= dcsol2(x,y,d);
            t0=toc;
            
            [b2,~ ,t2]= MMRN(x,y,initial,struct('d',d,'verbosity',0,'tolnorm',1e-7,'maxiter',1000));
            
            fun(k,:)=[DistCov(x*b0,y),DistCov(x*b2,y)];
            time(k,:)=[t0,t2];
            angle(k,:)=[subspace(b0,beta),subspace(b2,beta)];
        end       
        save(['ModelB_',num2str(i),'_',num2str(j)],'time','fun','angle');
    end
end

delete(gcp);
