function [Wmin,d,fmin] = permCORE(Yaux,X,morph,parameters)
%
% function [Wmin,d,f] = permLAD(Y,X,morph,parameters);
% 
% This function estimates the dimension of the reduced subspace that best
% describes the data under the CORE model using a permutation test.  
% 
% USAGE:
%  - outputs:
%    - Wmin: generating vectors for the central subespace of estimated
%    dimension.
%    - d: estimated dimension under LRT.
%    - f: value of the optimized function for dimension d. (perhaps this is useless)
%  - inputs: 
%    - Y: response vector;
%    - X: matrix of predictors;
%     morph: 'cont' for continuous responses or 'disc' for discrete
%     responses.
%     parameters (OPTIONAL): structure to set specific values of parameters for
%     computations. 
%           - parameters.nslices: number of slices for discretization of
%           continuous responses. % 5 slices are used by default.
%           - parameters.alpha: test level. Default is 0.05.
%           - parameters.npermute: number of permutations of the sample.
%           Default value is 500 permutations.
%           - parameters.sg: optional parameters for sg_min (see sg_min 
%           documentation for details).
%
%
%
% =========================================================================
%----checking type of response and slicing if needed.......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    if parameters.nslices==0,
        warning('MATLAB:slices','for continuous responses, a number of slices should be given. Five slices will be used');
        parameters.nslices = 5;
    end
    Y = slices(Yaux,parameters.nslices);
end

% ---- main process............................................................
h = parameters.nslices;
[n,p] = size(X);
data_parameters = setdatapars_v2(Y,X,h);

%--- get handle to objective function, derivative, dof ......................
Fhandle = F(@F4core,data_parameters);
dFhandle = dF(@dF4core,data_parameters);
% f0 = @(sigmag) (sum(sum(sigmag-sigmag)));

perm = get_perm_handle(Y,X,data_parameters,parameters);
d = perm(Fhandle,dFhandle);

if d>0
    [Wmin,fmin] = core(Y,X,d,morph,parameters);
else
    fmin = 0;
    Wmin = zeros(p);
end


