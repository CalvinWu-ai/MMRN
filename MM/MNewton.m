function dir=MNewton(Q,L,Y)
%  max 1/2*trace(Y'*Q*Y)+trace(Y'*L)
%  s.t. Y'*Y=I
% One Riemannian Newton's direction for the subproblem at the point Y.
% Q symmetry

% The code for the Algorithm 2 in the article
% We use Y to denote to the \gamma in the article

[p,d]=size(Y);
Yc=orthcompII(Y);
W1=Y'*Q*Y;
W2=Y'*Q*Yc;
W4=W2';
W5=Yc'*Q*Yc;
S=symm(W1+Y'*L);


I = speye(d,d);
I1= speye(p-d,p-d);
Td=Tvec(d);
Dd=veck_to_vec(d);

H11=kron(I,(W1-S))+kron((W1-S),I);
H11=(Dd'*H11*Dd)/4;

H12=kron(I,W2)-Td*kron(I,W2);
H12=1/4*Dd'*H12;
   
H21=kron(I,W4)*Dd;

H22=kron(I,W5)-kron(S,I1);

clear W1 W2 W3 W4 W5 Tp I I1

rgrad=projection(Y,Q*Y+L); % Riemannian gradient of the surrogate function at point Y
b1=-Y'*rgrad;
b1=b1(:);
b1=1/2*Dd'*b1;  % veck(U)=1/2*Dd'*vec(U)

b2=-Yc'*rgrad;
b2=b2(:);
b=[b1;b2];

H=[H11,H12;H21,H22];

clear H11 H12 H21 H22 b1 b2 rgrad 


x=H\b;
B=full(Dd*x(1:d*(d-1)/2));
B=reshape(B,d,d);

C=full(x((d*(d-1)/2+1):end));
C=reshape(C,p-d,d);

dir=Y*B+Yc*C; % Riemannian Newton's direction



end