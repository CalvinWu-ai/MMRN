function dcovXY = DistCov(x,y)

% DistCov calculates the distance covariance between the x vector
% and the y vector. 

% Input
% x: an n by p matrix
% y: an n by q matrix

% Output
% dcovXY: the distance covariance between the x vector and the
% y vector;

% References:
% [1] Gabor J. Szekely, Maria L. Rizzo and Nail K. Bakirov: Measuring and 
% Testing Dependence by Correlation of Distances
% Published by the Annals of Statistics
% 2007, Vol 35, No 6, 2769¨C2794


%------------------------------------------------initialize the calculation
sy=sum(y.^2,2);
B=real((bsxfun(@minus,sy',bsxfun(@minus,2*(y*y'),sy))).^(1/2)); % b_kl=||y_k-y_l||_2
B=B-mean(B,2)-mean(B,1)+mean(mean(B)); % normalization

sx=sum(x.^2,2);
A=real((bsxfun(@minus,sx',bsxfun(@minus,2*(x*x'),sx))).^(1/2)); % a_kl=||x_k-x_l||_2
%------------------------------------------------

dcovXY =mean(mean(A.*B));


end