
rng(1);
setpaths;

n=800;
p=20; 

P=[ones(1,p);repmat([1,-1],1,p/2)]'; % a p by 2 matrix
Phi=0.5.^abs((1:p)'-(1:p));

y=rand(n,1);
theta=2*pi*y; % theta~U[0,2*pi]
x=P*[sin(theta),cos(theta)]'+Phi*0.1*randn(p,n); % circle structure
x=x';

d=2;
initial=find_initial(x,y,d);

tic;
b0= dcsol2(x,y,d,initial);
t0=toc;

[b2,~ ,t2]= MMRN(x,y,initial,struct('d',d,'verbosity',0,'tolnorm',1e-7,'maxiter',1000));

coord1=x*b0;
coord2=x*b2; % a n by 2 matrix

subplot(121)
scatter(coord1(:,1),coord1(:,2),6,'filled');
axis([-2 2 -2 2]); axis square; box on; grid on;
xlabel('1st SDR component');ylabel('2nd SDR component');
title(['SQP (',num2str(t0,3) ' seconds)'])
subplot(122)
scatter(coord2(:,1),coord2(:,2),6,'filled');
axis([-2 2 -2 2]); axis square; box on; grid on;
xlabel('1st SDR component');ylabel('2nd SDR component');
title(['MMRN (',num2str(t2,3) ' seconds)'])






