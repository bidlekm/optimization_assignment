clear, clc;

data = load("Sharad.mat");
A = 10^-25;
a = data.a; dhdx = data.dhdx; g = data.g; rho = data.rho; x = data.x; H_obs = data.H_obs;

H = @(n) (-(2+n)/(2*A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^-1).^(1./(n+2));

sum1 = @(n) sum((H(n)- H_obs).^2);

fun = @(n) sqrt(sum1(n));

[nmin, fmin] = fminbnd(fun, 2, 3);
nmin
fmin

% Plot

i=0;
X = 2:1e-3:6;
F=zeros(1, length(X));
for n = X
    i = i + 1;
    F(i) = fun(n);
end

figure
plot(X, F)
xlim([1.5, 6])

%% Task 2
n = nmin;
H2 = @(A) (-(2+n)/(2*A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^-1).^(1./(n+2));

sum2 = @(A) sum((H2(A)- H_obs).^2);

fun2 = @(A) sqrt(sum2(A));

options = optimoptions('fminunc');
options.TolFun = 1e-28;
options.OptimalityTolerance = 1e-28;    % This completely sucks
options.StepTolerance = 1e-28;
options.MaxFunctionEvaluations = 5e3;
[Amin, fmin2] = fminunc(fun2, 2e-25, options);
Amin

i=0;
X2 = linspace(1e-28,1e-22,1e4);
F2=zeros(1, length(X2));
for n = X2
    i = i + 1;
    F2(i) = fun2(n);
end

figure
loglog(X2, F2)
ylim([5*1e2, 1e4])