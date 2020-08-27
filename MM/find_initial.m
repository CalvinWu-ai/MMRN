function b=find_initial(x,y,d)




[n,p]=size(x);
X=x-ones(n,1)*mean(x);
N=cov(X);
N1=sqrtm(N);
sX= x / N1;

F=@(beta) DistCov(sX*beta,y);

% Other Useful Initial Value Guess
%--------------------------------------------------------------------------
[~,W1,~] = SIR(y,sX,'cont',d,'nslices',6);
v1=orth(W1);
[~,W2,~] = DR(y,sX,'cont',d,'nslices',6);
v2=orth(W2);

b=v1;

if F(v2)>F(b)
    b=v2;
end

end



