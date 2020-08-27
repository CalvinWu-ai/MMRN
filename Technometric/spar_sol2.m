function [beta,st] = spar_sol2(G,N1,b0,lambda,itr)

% sparse solution using distance covariance. lam is the tuning parameter
% based on adaptive lasso. The weight fucntion is mywt(). itr is the number
% of iteration we want. The global variable wx is an intermediate variable used for removing
% redudant predictors. H is the matrix for LQA. The tolerance for
% convergence is set as 1e-3. 

global M;
global N;

N=N1;

A=sqrtm(N);
B=inv(A);
M=A*G*A;

tol=1e-3;

[p,d]=size(b0);
ww=b0;

st = ones(p,1);

if (lambda ==0)
    beta = b0;
    return;
end


v0=b0;

pw = mywt(ww);
H = diag(pw)*diag(1./sqrt(diag(v0*v0')));


for i=1:itr
   
   [V,D]=eig(G-0.5*lambda*B*H*B);
   [E,I] = sort(diag(D),'descend');
    V = V(:, I);
    b= B*V(:,1:d);
   dis = subspace(b0,b); % a function to calculate the largest principal angle for subspaces 
   v=b;
   [v,st,ix] = delp2(v,st); 
   
   if (sum(ix)>0)
    A=sqrtm(N);
    B=inv(A);
    G = B*M*B;
   ww(ix,:)=[];
   b(ix,:)=[];
   pw = mywt(ww); 
   end
   
   if (sum(ix)==0)&&(dis <= tol)
       break;
   end

   H = diag(pw)*diag(1./sqrt(diag(v*v')));
   if (size(H,1)<=d) 
       break;
   end
   b0=b;
end


te = 0;
for i=1:p
    if st(i)==1
        te = te+1;
    beta(i,1:d) = v(te,1:d);
    elseif st(i)==0
        beta(i,1:d) = zeros(1,d);
    end
end





  



    
       