function [nx,index]=seq_dc(x,y,m,H)

[n,p]=size(x);

p1=20;
itr=200;
index=1:p;

x1=x(:,1:p1);
x2=x(:,p1+1:p);
wy=[x2,y];
ncol=p;

while (ncol > m)

        [b0,d,G,N] = PRSIR(wy,x1,H);
        if (d==0)
            x(:,1:p1)=[];
            ncol=ncol-p1;
            index(1:p1)=[];
            if (ncol<p1)
                 nx=x;
                 break
             end
            x1=x(:,1:p1);
            x2=x(:,p1+1:ncol);
            wy=[x2,y];
        else
            [beta,st] = bic_sol2(n,G,N,b0,itr);
            if (sum(st)==p1)
            nx=x;
            break
            end

            k=1;
            tx=[];
            for i=1:p1
                if (st(i)==0)
                    tx(k)=i;
                    k=k+1;
                end
            end
            

            
            index(tx)=[];
            x(:,tx)=[];
            [temp,ncol]=size(x);
            
            
           if (ncol <m)
           nx=x;
           break
           end
            
            x1=x(:,1:p1);
            x2=x(:,p1+1:ncol);
            wy=[x2,y];
        end
        
end

nx=x;
