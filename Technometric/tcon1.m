function [c,ceq]  = tcon1(b)

% Constraint function using the covariance matrix

global N;

d=size(b,2);

c = max(abs(b))-5;     % Compute nonlinear inequalities at x.
ceq = b'*N*b-eye(d); % Compute nonlinear equalities at x.