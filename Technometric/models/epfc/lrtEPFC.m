function [Wmin,d,f] = lrtEPFC(Yaux,X,morph,parameters)
%
% [Wmin,d,f] = lrtEPFC(Y,X,morph,parameters);
% 
% This function estimates the dimension of the central subspace  under 
% the EPFC model using a likelihood-ratio test.
% USAGE:
%  - outputs:
%    - Wmin: generating vectors for the central subespace of estimated
%    dimension.
%    - d: estimated dimension under LRT.
%    - f: value of the optimized function for dimension d. (perhaps this is useless)
%  - inputs: 
%    - Y: response vector;
%    - X: matrix of predictors;
%    - morph: 'cont' for continuous responses or 'disc' for discrete
%     responses.
%    - parameters (OPTIONAL): structure to set specific values of parameters for
%     computations. 
%           - parameters.alpha: test level. Default is 0.05.
%           - parameters.fy: basis for regression of centered predictors.
%           By default, this is a polynomial basis of degree 3.
%
%
%
% -------------------------- REQUIREMENTS ---------------------------------
% This function requires the Statistic Toolbox or a custom function to 
% compute the CDF for a chi2 distribution.
% =========================================================================

% ====================================================================

%----checking type of response ......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    Y = Yaux;
    parameters.nslices = length(Y);
end        

%----read size of sample 
[n,p] = size(X);

%--- get sample statistics ................................................
data_parameters = setdatapars_v2(Y,X,parameters.nslices);
if strcmpi(morph,'cont')
    [SIGMAfit,r] = get_fitted_cov(Y,X,parameters.fy);
else
    [SIGMAfit,r] = get_average_cov(X,data_parameters);
end
SIGMA = data_parameters.sigmag;
SIGMAres = SIGMA - SIGMAfit;
data_parameters.Afit = SIGMAfit;
data_parameters.B = SIGMAres;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4epfc,data_parameters);
dFhandle = dF(@dF4epfc,data_parameters);
dof = @(u)(r*u + p*(p+3)/2);
f0 = n*p/2 + n*p/2*log(2*pi);

%--- slices response Y to get initial estimates using SIR, etc
if strcmpi(morph,'cont')
    haux = 5;
    Ysliced = slices(Y,haux);
    aux_datapars = setdatapars_v2(Ysliced,X,haux);
    auxpars = parameters; auxpars.nslices=haux;
else
    Ysliced = Y;
    aux_datapars = data_parameters;
    auxpars = parameters;
end


lrt = lrtloop_handle_vf(Y,X,data_parameters,parameters,Ysliced,aux_datapars,auxpars);
d = lrt(Fhandle,dFhandle,f0,dof);

if d>0
    [Wmin,f] = epfc(Y,X,d,morph,parameters);
else
    f=f0;
    Wmin=zeros(p);
end