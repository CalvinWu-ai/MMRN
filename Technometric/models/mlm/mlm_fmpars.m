function [beta,E,S,Sfit,Sres] = mlm_fmpars(Y,X);
% [beta,S,Sfit,Sres] = mlm_fmpars(Y,X)
% Returns the parameter estimates from the full multivariate linear model.
% beta is the coefficient matrix
% S is the marginal covariance matrix of Y
% Sfit is the covariance matrix of the fitted vectors
% Sres is the covariance matrix of the residual vectors
% ====================================================
n = size(X,1);
p = size(X,2);
F = mlm_center(X);
U = mlm_center(Y);
beta = U'*F*inv(F'*F);
E=U-F*beta';
Sfit = beta*F'*U/n;
S = U'*U/n;
Sres = S - Sfit;
