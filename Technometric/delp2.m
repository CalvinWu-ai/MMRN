function [nb,nst,ix] = delp2(b,st)

% the function used for removing a predictor when its norm of corresponding
% coefficient is small enough. The tolerance level is set as 1e-3. The parameter st is to monitor the status of
% predictors. 

global M;
global N;
tol=1e-4;

p=length(st);

t = 0;
s = 0;
nst=st;

ix=0;

for i=1:p
    
   if st(i) == 0  
       continue
   else
       t = t +1;
       if norm(b(t,:)) < tol
           s = s+1;
           nst(i)=0;
           ix(s)=t;
       else
           continue
       end
   end
end

if sum(ix)>0
b(ix,:)=[];
M(:,ix)=[];
M(ix,:)=[];
N(:,ix)=[];
N(ix,:)=[];
end

nb=b;
