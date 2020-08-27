function [dcovXY,dcorrXY,dvarXX, dvarYY] = DistCorrVec(x,y)

% DistCorrVec calculates the distance correlation between the x vector
% and the y vector. 

% Input
% x: an n by p matrix
% y: an n by q matrix

% Output
% dcorrXY: the distance correlation between the x vector and the
% y vector;
% dcovXY: the distance covariance between the x vector and the
% y vector;
% dvarXX: the distance variance of the x vector;
% dvarYY: the distance variance of the y vector.

% References:
% [1] Gabor J. Szekely, Maria L. Rizzo and Nail K. Bakirov: Measuring and 
% Testing Dependence by Correlation of Distances
% Published by the Annals of Statistics
% 2007, Vol 35, No 6, 2769¨C2794

% A Toy Example:
% n = 500; p = 1; q = 1;
% x = normrnd(0,1,n,p); y = x(:,1:q);
% dcorrXY = DistCorrVec(x,y)


[n,p] = size(x); II = ones(n,1); 
sxy1 = zeros(n,1); sxy2 = zeros(n,1); sxy3 = zeros(n,1); 
sxx1 = zeros(n,1); syy1 = zeros(n,1); 
%------------------------------------------------initialize the calculation
for ii = 1:n
    XX1 = sqrt(sum((x - II*x(ii,:)).^2,2));
    YY1 = sqrt(sum((y - II*y(ii,:)).^2,2));
    sxy1(ii) = mean(XX1 .* YY1);
    sxy2(ii) = mean(XX1);
    sxy3(ii) = mean(YY1);

    XX2 = XX1.^2; sxx1(ii) = mean(XX2);
    YY2 = YY1.^2; syy1(ii) = mean(YY2);
end
%---------------------------------------utilize moment estimator throughout

SXY1 = mean(sxy1);
SXY2 = mean(sxy2)*mean(sxy3);
SXY3 = mean(sxy2 .* sxy3); 
%-----------------------------------------calculate the distance covariance

SXX1 = mean(sxx1);
SXX2 = (mean(sxy2)).^2;
SXX3 = mean(sxy2.^2);
%-------------------------------------------calculate the distance variance

SYY1 = mean(syy1);
SYY2 = (mean(sxy3)).^2;
SYY3 = mean(sxy3.^2);
%-------------------------------------------calculate the distance variance

dcovXY = sqrt(SXY1 + SXY2 - 2 * SXY3);
dvarXX = sqrt(SXX1 + SXX2 - 2 * SXX3);
dvarYY = sqrt(SYY1 + SYY2 - 2 * SYY3);
dcorrXY = dcovXY ./ sqrt(dvarXX * dvarYY);
%----------------------------------------calculate the distance correlation
