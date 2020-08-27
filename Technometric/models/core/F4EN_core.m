function f = F4EN_core(W,FParameters)
% Objective function (minus the likelihood) for the CORE model.
% USAGE: 
% - W is the projection matrix onto the reduced subspace.
% Notice that global FParameters are supposed to be already set.
% ==============================================================
sigma = FParameters.sigma;
sigmag = FParameters.sigmag;
n = FParameters.n;
% ---define some convenience variables
h = n/sum(n);
a = zeros(length(h),1);
sigmatmp = zeros(size(sigma,2),size(sigma,3));
for i=1:length(h),
    sigmatmp(:,:) = sigma(i,:,:);
    a(i) = logdet(W'*sigmatmp*W);
end
% ---Likelihood function for LAD model
f = sum(n)/2 * (logdet(W'*inv(sigmag)*W) + h*a);
