function Gamma0 = mlm_Gamma0(Gamma);
% Gamma0 = mlm_Gamma0(Gamma)
% Returns Gamma0 given Gamma
% ===========================
r = size(Gamma,1);
u = size(Gamma,2);
QGamma = eye(r) - Gamma*inv(Gamma'*Gamma)*Gamma';
%[Gamma0,D] = firsteigs(QGamma,r-u);
Gamma0 = firsteigs(QGamma,r-u);
