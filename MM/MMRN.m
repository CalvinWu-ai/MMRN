function [beta,object,time] = MMRN(x,y,gamma0,option)
% MM Riemannian Newton (MMRN) Algorithm for solving DCOV based SDR model.
% max V^2(\beta'X,Y)
% s.t. \beta'cov(X)*\beta=I_d

%
% function [beta] = MMRN(x,y)
% function [beta] = MMRN(x,y,[],option)
% function [beta] = MMRN(x,y,gamma0,option)
%

% Input
% x: an n by p matrix
% y: an n by q matrix
% gamma0: a p by d orthogonal initial value matrix
% option: The options structure is used to overwrite the default values. All
%         options have a default value and are hence optional.
%         tolnorm (1e-7): The algorithm terminates if consecutive function ratio less. 
%         maxiter (1000): The algorithm terminates if maxiter iterations were executed.
%         epsilon (1e-10); The perturbed constant.
%         d(1): The dimension number of central subspace 
%         verbosity (2): Integer number used to tune the amount of output the algorithm
%                        generates during execution (mostly as text in the command window).
%                        The higher, the more output. 0 means silent. 2 and above includes a
%                        display of the options structure at the beginning of the execution.

% Output
% beta: the solution
% object: store the perturbed objective value at each iterate
% time: CPU time
%


[n,p]=size(x);


%------------------------------------------------defaults for option
maxiter=1000;
tolnorm=1e-7;
d=1;
epsilon=1e-10; % Perturbed constant
verbosity=2; % verbosity>=2 display every iteration message, otherwise display the final convergent message
%-------------------------------------------------



% User's option
%------------------------------------------------
if nargin > 3  %  then paramstruct is an argument
  if isfield(option,'maxiter')    
    maxiter = option.maxiter; 
  end 
  if isfield(option,'tolnorm')   
    tolnorm = option.tolnorm; 
  end  
  if isfield(option,'d')   
    d = option.d; 
  end
  if isfield(option,'epsilon')   
    epsilon = option.epsilon; 
  end
  if isfield(option,'verbosity')   
    verbosity = option.verbosity; 
  end
end
%---------------------------------------------------

Kermat=@(x) Kermatrix(x);
Newton=@(A,B,W0)  MNewton(A,B,W0); % Riemannian Newton's direction

B=Kermat(y);
B=B-mean(B,2)-mean(B,1)+mean(mean(B)); 


N=cov(x);
N1=sqrtm(N);
z=x/N1;   % Standardized x
Index=Kermat(z);
Index=(Index~=0);

F=@(gamma) PDistCov(z*gamma,y,epsilon); % Perturbed objective function

% Create a random starting point if no starting point is provided
if ~exist('gamma0', 'var')|| isempty(gamma0)
    gamma = orth(randn(p,d)); % random generating a p*d matrix   
else
    gamma = gamma0;
end



tic;

iter=0;
F_cur=F(gamma);


if verbosity >= 2
    fprintf('   iter   sub_iter    cost val \t  \t \t      tolnorm  \n');
end

% main code
error=1;
while true
    if iter == maxiter
        if verbosity>0
            disp('MMRN terminates: Achieved maximum iteration.');
        end
        break;
    end
    
    if  error < tolnorm
        if verbosity>0
            fprintf('MMRN terminates: converged iteration:%4d\n', iter);
        end
        break;
    end
   
    
    % Coefficients of Quadratic Term in Surrogate Function
    xgamma=Kermat(z*gamma);
    xgamma(Index)=1./(xgamma(Index)+epsilon);
    xgamma=xgamma.*Index;
    C=(B.*(B<0)).*xgamma;
    Q=z'*(2*(diag(sum(C,2))-C)/(n^2))*z;
    
    % Coefficients of Linear Term in Surrogate Function
    D=(B.*(B>0)).*xgamma;
    L=z'*(2*(diag(sum(D,2))-D)/(n^2))*z*gamma;
   
    
    
   %% Riemannian Newton Method for Solving the Subproblem
    
    gamma_pre=gamma;
    dir=Newton(Q,L,gamma_pre); 
    
    step=1;
    [gamma,~]=qr(gamma_pre+step*dir,0); % qr Retr operation
    F_trial=F(gamma);
    normDsquare=trace(dir'*dir);
    
   %% linesearch
    sub_iter=0;
    while F_trial <= F_cur+(1e-20)*step*normDsquare
        step=0.5*step;
        sub_iter=sub_iter+1;
        if step < 1e-10
            break;
        end
        [gamma,~]=qr(gamma_pre+step*dir,0);
        F_trial=F(gamma);
    end
    %%
      
    iter=iter+1;
    F_trial=F(gamma);
    object(iter)=F_trial;    
    
    
    % Display iteration information
    error=abs((F_trial-F_cur)/F_cur);    
    
    if verbosity >= 2
        fprintf('%5d  \t %5d \t  %.8e  \t      %.8e   \n', ...
                iter, sub_iter, F_trial, error);
    end
    
    
    F_cur=F_trial;

end

time=toc;
if verbosity >= 1
    fprintf('Total time is %f [s] \n', ...
        time);
end
beta=N1\gamma;

end



function K=Kermatrix(x)

% Input
% x: an n by p matrix
% 
% Output
% K: an n by n matrix with B_{ij}=||x_i-x_j||_2


sx=sum(x.^2,2);
K=real((bsxfun(@minus,sx',bsxfun(@minus,2*(x*x'),sx))).^(1/2)); % K_ij=||x_i-x_j||_2

end


function dcovXY = PDistCov(x,y,ep)

% Input
% x: an n by p matrix
% y: an n by q matrix

% Output
% dcovXY: the perturbed distance covariance between the x vector and the
% y vector;


%------------------------------------------------initialize the calculation
sy=sum(y.^2,2);
B=sqrt(real((bsxfun(@minus,sy',bsxfun(@minus,2*(y*y'),sy))))); % b_kl=||y_k-y_l||
B=B-mean(B,2)-mean(B,1)+mean(mean(B)); % normalization

sx=sum(x.^2,2);
A=sqrt(real((bsxfun(@minus,sx',bsxfun(@minus,2*(x*x'),sx))))); % a_kl=||x_k-x_l||
A=A-ep.*log(1+A./ep); % Perturbed Objective Function
%------------------------------------------------

dcovXY =mean(mean(A.*B));


end




