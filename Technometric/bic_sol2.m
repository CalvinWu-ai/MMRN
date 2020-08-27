function [beta,st] = bic_sol2(n,G,N,b0,itr)

[p,d]=size(b0);

A=sqrtm(N);
B=inv(A);
M=A*G*A;


for i=1:31
     [beta,st] =spar_sol2(G,N,b0,i/100-0.01,itr);
     f = -trace(beta'*M*beta) + 2*d*(sum(st)-d)/n;
    stos(i) = f;
end

[f,ind] = min(stos);
las = ind/100-0.01;

for k = 1:19
    

    [beta,st] =spar_sol2(G,N,b0,las-0.01+k/1000,itr);
    f = -trace(beta'*M*beta) + 2*d*(sum(st)-d)/n;
    dstos(k) = f;
end

[f,dind] = min(dstos);
dlas = las-0.01+dind/1000;

[beta,st] = spar_sol2(G,N,b0,dlas,itr);


