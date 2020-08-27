function B = dMAVE(x, y, m0);

% Input : x (n x p); y (n x 1); m0 (integer)
% Output: B (p x m0) 


[n,p] = size(x);
onen = ones(n,1);

mm = mean(x);
x = x - repmat(mm,n,1);
ss = inv(x'*x/n)^0.5;
x = x*ss;

m = p;
B = eye(p);
Ba = B;
BI = Ba;
B = B(:,1:m);
noip0 = 0;
noip1 = 1;
iterstop = 0;
Btmp = B;
rige = std(y)*mean(std(x,[],1));
y = (y-mean(y))/std(y);
x = x./repmat(std(x,[],1),n,1);

yc = y;
x0 = x;
xc = x;

if 1 == 1
[a, yt] = sort(y);
[a, yt] = sort(yt);
ycc = yt/max(yt);
for i = 1:size(x,2);
    a = x(:,i);   
    [b, a] = sort(a);
    [b, a] = sort(a);
    xc(:,i) = a/std(a);
end
end

niter = floor(p*3/2);
%ch = (sqrt(p)/n^(1/(p+4))*n^(2/(m0+4)))^(1/niter);
for iter = 1:niter;
    x = xc;
    if iter >= p;
        x = x0;
    end
  	adj = p^(1/iter)*n^(1/(m0+4)-1/(m+4));
  	hy = std(yc)/n^(1/5)*p^(1/(iter+1)); %*adj;
	ky = repmat(yc, 1, n);
    ky = [ycc/sqrt(2*pi)/hy*n^(1/iter) exp(-(ky-ky').^2/(2*hy*hy))/sqrt(2*pi)/hy];
    n1 = size(ky,2);
    h = p^(1/(iter+1))*mean(std(x*Ba,[],1))/n^(1/(m+4));
    h2 = 2*h*h*adj;
   
   ABI = zeros(m, n*n); 
   for iiter = 1:max(1, (m < p)*floor(m/2));
	dd = zeros(m*p, m*p);
   	dc = zeros(m*p,1);
   	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + sum(xij.^2,2)/1.5^iter;
      	ker = exp(-sxij/h2);  
        rker = repmat(ker, 1, p+1);
   		onexi = [xij*B onen];
   		xk = (onexi.*rker(:, 1:m+1))';
   		abi = inv(xk*onexi+eye(m+1)/n)*xk*ky;
      
      	kxij = (xij.*rker(:,1:p))';
      	kxijy = kxij*(ky - repmat(abi(m+1,:),n,1));
      	ddx = kxij*xij;
      	for k1 = 1:m
         	ka = (k1-1)*p+1:k1*p;
         	dc(ka) = dc(ka) + kxijy*abi(k1,:)';
      	end
        tmp = abi(1:m, :)*abi(1:m, :)';
        dd = dd + kron(tmp, ddx);
      
        ABI(:, (j-1)*n1+1:j*n1) = abi(1:m,:); 
   	end
    
	   B = pinv(dd+ rige*eye(length(dc))/n)*dc;
       
   	   B0 = reshape(B, p, m)*ABI;
       [B, R] = eig(B0*B0');
       B = B(:,p-m+1:p);
       Ba = B;
       if (max(svd(B*B'-BI*BI')) < 0.001);
           break
       end
       BI = B;
       
   end

   mb = m;
   ma = max(m-1,m0);
   m = ma; 
   B = B(:,mb-m+1:mb);
   Ba = B;
   if (max(svd(B*B'-BI*BI')) < 0.001)*(iter>p+3);
       break
   end
   BI = Ba;   
   
end

cv = 0;
kye = ky;
	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + 0*sum(xij.^2,2)/1.5^iter;
        ker = exp(-sxij/h2);  
        ker(j) = 0;
        if mean(ker) > 1/n/n;
           kye(j,:) = ker'*ky/sum(ker);
        end      
  end   
cv = (ky-kye).^2;
cv = sum(sum(cv))/sum(sum(cv>0));

B = ss*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end
