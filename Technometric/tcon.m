function [c,ceq]  = tcon(x)

% Usual constraint function

d=size(x,2);

c = [];     % Compute nonlinear inequalities at x.
ceq = x'*x-eye(d); % Compute nonlinear equalities at x.