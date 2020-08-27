function [beta,st] = spar_sol(x,y,b0,lam,itr)

% sparse solution using distance covariance. lam is the tuning parameter
% based on adaptive lasso. The weight fucntion is mywt(). itr is the number
% of iteration we want. The global variable wx is an intermediate variable used for removing
% redudant predictors. H is the matrix for LQA. The tolerance for
% convergence is set as 1e-3. 

global wx;
global H; 
global lambda
global N2;


wx=x;

tol=1*1e-3;

[p,d]=size(b0);

ww=b0;

st = ones(p,1);

lambda = lam;

N=cov(wx);

N1=sqrtm(N);

N2=inv(N1);

if (lambda ==0)
    beta = b0;
    return;
end

options = optimset('Display','off','Algorithm','sqp','MaxFunEvals',1*1e6,'MaxIter',1*1e6,'TolFun',1*1e-4,'TolX',1*1e-4,'LargeScale', 'off');



pw = mywt(ww);


[b0,st,ix] = delp(b0,st); 
 
   if (sum(st)<d)
       disp('something wrong');
       return;
   end
   
   
   if (sum(ix)>0)
   ww(ix,:)=[];
   pw = mywt(ww); 
   N=cov(wx);
   N1=sqrtm(N);
   N2=inv(N1);
   end
   
v0=N1*b0;

H = diag(pw)*diag(1./sqrt(diag(b0*b0')));

for i=1:itr
    
   
   [v,fv,exitflag] = fmincon(@(c)rdcreg(c,wx,y),v0,[],[],[],[],[],[],@tcon,options);
   
   if (exitflag==-1)
       beta=zeros(p,1);
       st=zeros(p,1);
       return;
   end
   
   
   b= N2*v;
   
   dis = subspace(b0,b); % a function to calculate the largest principal angle for subspaces 
   
   [nb,nst,ix] = delp(b,st); 
   
   if (sum(nst)<d)
       break;
   end
   
   b=nb;
   st=nst; 
   
   if (sum(ix)>0)
   ww(ix,:)=[];
   pw = mywt(ww); 
   N=cov(wx);
   N1=sqrtm(N);
   N2=inv(N1);
   end
   

   
   if (sum(st)==d)
       ix=0;
   end
   
   
   if (sum(ix)==0)&&(dis <= tol)
       break;
   end

   
   H = diag(pw)*diag(1./sqrt(diag(b*b')));

   b0=b;
   v0=N1*b0;
end


te = 0;
for i=1:p
    if st(i)==1
        te = te+1;
    beta(i,1:d) = b(te,1:d);
    elseif st(i)==0
        beta(i,1:d) = zeros(1,d);
    end
end

wx=x;

return;




  



    
       