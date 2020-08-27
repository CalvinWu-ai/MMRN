function f= rdcreg(v,wx,y)

% the appoximated objective function using LQA. 

global H; 
global lambda
global N2;

g=DistCorrVec(wx*N2*v,y);

f=-g + 0.5*lambda*trace(v'*N2*H*N2*v);

return;