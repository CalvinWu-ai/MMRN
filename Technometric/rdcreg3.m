function f= rdcreg3(b,wx,y)

% the appoximated objective function using LQA. 

global H; 
global lambda


g=DistCorrVec(wx*b,y);

f=-g + 0.5*lambda*trace(b'*H*b);

return;