function [beta,st,time,stos] = PMMRNBic(x,y,b0,rang)
% DCOV based SVS tunning parameter using BIC criterion
%
% function [beta,st,time,stos] = PMMRNBic(x,y,b0,rang)
%
% Input
% x: an n by p matrix
% y: an n by q matrix
% b0: the initial value, usually the estimate without penalty
% rang: the range for searching the solution. Suggestion area from 0.5 to 2;
%
% Output
% beta: the solution
% st: the selection status of variables
% time: the total running time
% stos: BIC value at each tunning parameter

% st:  the selection status of variables
tic;
[n,~]=size(x);

[~,d]=size(b0);

IR=floor(rang*100)+1;
theta=sqrt(sum(b0.^2,2)).^(-0.5);

for i=1:IR % the searching range 
    [beta,st] = PMMRN(x,y,b0,struct('theta',theta,'d',d,'verbosity',0,'Lambda',i/100-0.01));
    g=DistCov(x*beta,y);
    f = -log(g) + log(n)*d*(sum(st)-d)/n;   % BIC   
    %f=-log(g)+2*d*(sum(st)-d)/n; % AIC
    stos(i) = f;
end

[~,ind] = min(stos);
las = ind/100-0.01;

for k = 1:19 % the searching sub-range 
    
    [beta,st] = PMMRN(x,y,b0,struct('theta',theta,'d',d,'verbosity',0,'Lambda',las-0.01+k/1000));
    g=DistCov(x*beta,y);
    f = -log(g) + log(n)*d*(sum(st)-d)/n;   % BIC  
    %f=-log(g)+2*d*(sum(st)-d)/n; % AIC
    dstos(k) = f;  
    
end

[~,dind] = min(dstos);
dlas = las-0.01+dind/1000;



[beta,st] = PMMRN(x,y,b0,struct('theta',theta,'d',d,'verbosity',0,'Lambda',dlas));

time=toc;
end

