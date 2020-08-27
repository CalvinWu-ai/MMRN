data = load('boston_trim.txt');

y = data(:,1);
x = data(:,2:14);

d=2;
rang=0.5;
setpaths;

b0= dcsol2(x,y,d);
[beta,st] = bic_sol(x,y,b0,rang,500);