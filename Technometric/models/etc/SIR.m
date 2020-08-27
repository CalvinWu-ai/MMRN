function [WX,W,auxmtx] = SIR(Yaux,X,morph,dim,varargin)
% function [WX,W] = SIR(Y,X,morph,u,varargin)
% 
% This function implements the sliced inverse regression procedure for 
% sufficient dimensionality reduction. See references for details.
%
% USAGE:
%   - outputs:
%     WX: projection of the predictors onto the dimension reduction subspace.
%     W: generating vectors of the dimension reduction subspace.
%   - inputs:
%     Y: response vector.
%     X: predictors matrix.
%     morph: with value 'cont', specifies that the response Y is continuous 
%     (in which case it is a regression problem) while with value 'disc' it
%     specifies a discrete response (and a classification problem).
%     u: dimension for the reduced subspace.
%     varargin: optional arguments. They must be given as a string-number pair in which
%     the string specifies de optional parameter and the following input sets its value.
%     Available options are limited to:
%     - 'nslices': to set the number of slices to be used to discretize continuous 
%       responses.
%     - 'setmtx': boolean flag used to specify if the computation will rely on previous
%       auxiliary results. This is useful when several sufficient dimensionality 
%       reduction methods are tested. These methods often share some procedures
%       that can be performed before computation to speed up the process. 
%       As a flag, allowed values are:
%       - true or 1: to use auxiliary results stored in a global variable by 
%                    using function SETAUX.
%       - false or 0: to perform all computations.

% =========================================================================


%----checking required arguments...........................................
if nargin < 4,
    error('Not enough input arguments. Type >>help ldr for details');
end
% .........................................................................

%----checking data consistency.............................................
if size(Yaux,1)~=size(X,1),
    error('The number of rows in Y must equate the number of rows in X');
end
if ~strcmpi(morph,'cont') && ~strcmpi(morph,'disc'),
    error('unknown type of response. Valid options are CONT or DISC...');
end
if ~ischar(dim) && ~isinZ(dim),
    error('Natural value expected to specify reduced subspace dimension.');
end
if ~isreal(Yaux),
    error('Response vector must be numeric');
end

%----reading optional input arguments and saving parameters................
% if nargin > 4,
    parameters = read_input_nldr(varargin{:});
% else
%     parameters.auxhandle
% end

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


sigmag = get_cov(X);
p = size(sigmag,2);
[V,D] = eig(sigmag);
if isempty(parameters.auxhandle),
    parameters.auxhandle = setaux_v2(X,Y,parameters.nslices,p,V,D);
end
auxmtx = parameters.auxhandle;
[SIRdir,SIR] = getSIRv2(dim,parameters.auxhandle);

%-----Write results------------------------
W = SIR;
WX = X*W;
    
    
