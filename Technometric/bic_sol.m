function [beta,st] = bic_sol(x,y,b0,rang,itr)

% rang is the range for searching the solution. Suggestion area from 0.5 to 2;
% b0 is the initial value, usually the estimate without penalty
% st is the selection status of variables

[n,p]=size(x);

[p,d]=size(b0);

IR=floor(rang*100)+1;

for i=1:IR % the searching range 
    [beta,st] = spar_sol(x,y,b0,i/100-0.01,itr);
    g=DistCorrVec(x*beta,y);
    f = -2*log(g) + log(n)*d*(sum(st)-d)/n;
    stos(i) = f;
end

[f,ind] = min(stos);
las = ind/100-0.01;


for k = 1:19 % the searching sub-range
    
    [beta,st] = spar_sol(x,y,b0,las-0.01+k/1000,itr);
    
    g=DistCorrVec(x*beta,y); 
    
    f = -2*log(g) + log(n)*d*(sum(st)-d)/n;
    dstos(k) = f;
end

[f,dind] = min(dstos);
dlas = las-0.01+dind/1000;


[beta,st] = spar_sol(x,y,b0,dlas,itr);

