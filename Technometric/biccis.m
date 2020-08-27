function [f,beta,st] = biccis(lam)

global lambda;
global n;
global d;
global p;
global M;
global N;  

M0= M;
N0= N;

if lam < 0
    f = inf;
    return
end

lambda = lam;

[A,B] = mysqrt(N);
G = B*M*B;
[a2,b2]= eig(G);
[a1,b1] = sort(diag(b2));

if d == 2
x0 = B*[a2(:,b1(p)),a2(:,b1(p-1))];

elseif d==1
    x0 = B*a2(:,b1(p));
end

[beta,st] = ecis(x0);

M = M0;
N = N0;



f = -trace(beta'*M*beta) + d*(sum(st)-d)*log(n)/n;


