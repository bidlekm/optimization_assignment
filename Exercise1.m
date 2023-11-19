clear, clc;
data = load("D:\Uppsala\Optimization\Sharad.mat");
A = 10^-25;
n = balance_function(data.a,data.dhdx,data.g,data.rho,data.x,A,data.H_obs)