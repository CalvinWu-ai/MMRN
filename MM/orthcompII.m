function YCOMP=orthcompII(Y)
% Compute the orthogonal complement of Y

[n,p]=size(Y);
[Q,~] = qr(Y);
YCOMP=Q(:,(p+1):n);
end