function [beta,SI] = ESE(x,y,d,rang,itr)
% This is the main function
% d is the dimension of the central subspace 
% itr is the number of interation in sparse estimating
% SI is the index of selected important predictors
% beta is the correspoding coefficients
% rang is the range for searching the solution. Suggestion area from 0.5 to 2;

m=30;% the upper size of last step when stopping the sequential process 
H=6; % number of slices for the projection resampling SIR
setpaths;

[n,p]=size(x); % p is the number of predictors, n is the number of obersvation
[n,r]=size(y); % r is the dimension of the response

if (r==1)  % to justify if the response is univariate or not
    if (n>p)  
    b0= dcsol2(x,y,d); % initial estimate of the central subspace 
    [beta,st] = bic_sol(x,y,b0,rang,itr);
    L=logical(st);
    K=1:p;
    SI=K(L);
        else
        [nx1,index1] = dc_ordering(y,x); % sorting variables from low to high distance correlation
        [nx2,index2]=seq_dc1(nx1,y,m,H); % sequential sparse dimension reduction
        idx=index1(index2);
        b0= dcsol2(nx2,y,d);
        [beta,st] = aic_sol(nx2,y,b0,itr);% final step of sparse estimation
        L=logical(st);
        SI=idx(L);
    end
else 
    b0 = PRSIR1(y,x,d,6);% initial estimate of the central subspace using projection resampling for mulitvariate response
    % at this moment, we require n>p when r>1. 
    [beta,st] = bic_sol(x,y,b0,itr);  
    L=logical(st);
    K=1:p;
    SI=K(L);
    
end

end

