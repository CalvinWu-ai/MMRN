function df = dF4EN_core(W,FParameters)
%	Derivative of F (minus the likelihood) on the dimension of the 
%   reduced subspace for the CORE model.
%   USAGE: 
%     - W is the projection matrix onto the reduced subspace, as
%       computed with the FNFP function.
%   Notice that global FParameters are supposed to be already set.
% ==============================================================
    sigma = FParameters.sigma;
    sigmag = FParameters.sigmag;
    n = FParameters.n;
    f = n/sum(n);
    a = zeros(size(W,1), size(W,2), length(f));
    firsta = zeros(size(sigma,2),size(sigma,3));
    for i=1:length(f)
        firsta(:,:) = sigma(i,:,:);
        a(:,:,i) = sum(n)*f(i)*firsta*W*inv(W'*firsta*W);
    end
    first = zeros(size(W));
    for j=1:length(f)
        first = first+a(:,:,j);
    end
    second=sum(n) * inv(sigmag) * W * inv(W' * inv(sigmag) * W);
    % ---Derivative of likelihood function for the LAD model
    df = (first+second);
