clear; clc;

data = load("Sharad.mat");
A = 10^-23;
a = data.a; dhdx = data.dhdx; g = data.g; rho = data.rho; x = data.x; H_obs = data.H_obs;

H = @(n, A) (-(2+n)/(2*10.^A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^(-1)).^(1./(n+2));

sum1 = @(args) H(args(1), args(2))- H_obs;

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');

sol = lsqnonlin(sum1, [2, -24],[],[], options)