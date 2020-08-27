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
        beta(:,1)=[1,0.5,1,0,0,0,zeros(1,p-6)]';
        d=1;
                
        time=zeros(100,3);
        angle=zeros(100,3);
        fun=zeros(100,3);
        parfor k=1:100
            if j==1
                x=randn(n,p);
            elseif j==2
                x=betarnd(1.5,1,[n p]); 
                x=x.*2-1;
            else
                x=zeros(n,p);
                x(:,1:5)=poissrnd(1,[n 5]);
                x(:,6)=binornd(10,0.3,[n,1]);
                x(:,7:p)=poissrnd(1,[n p-6]);
            end
            temp=x*beta;
            y=exp(temp).*randn(n,1);
            initial=find_initial(x,y,d);

            tic;
            b0= dcsol2(x,y,d);
            t0=toc;
            
            [b2,~ ,t2]= MMRN(x,y,initial,struct('d',d,'verbosity',0,'tolnorm',1e-7,'maxiter',1000));
           
            fun(k,:)=[DistCov(x*b0,y),DistCov(x*b2,y)];
            time(k,:)=[t0,t2];
            angle(k,:)=[subspace(b0,beta),subspace(b2,beta)];
        end       
        save(['ModelC_',num2str(i),'_',num2str(j)],'time','fun','angle');
    end
end

delete(gcp);
