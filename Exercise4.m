clear, clc;

data = load("Sharad.mat");
A = 10^-23;
a = data.a; dhdx = data.dhdx; g = data.g; rho = data.rho; x = data.x; H_obs = data.H_obs;

H = @(n, A) (-(2+n)/(2*A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^-1).^(1./(n+2));

sum1 = @(n, A) sum((H(n,A)- H_obs).^2);

fun = @(args) sqrt(sum1(args(1),args(2)));
initialdata = [1.1, 1e-24];



options = optimoptions('fminunc');
options.TolFun = 1e-35;
options.OptimalityTolerance = 1e-35;    % This completely sucks
options.StepTolerance = 1e-35;
options.MaxFunctionEvaluations = 5e4;
[minParams, fval] = fminunc(fun, initialdata, options);
minParams
fval







n_values = linspace(2, 3, 100);
A_values = linspace(1e-22, 1e-26, 100);

[NN, AA] = meshgrid(n_values, A_values);
Z = zeros(size(NN));

for i = 1:numel(NN)
    Z(i) = fun([NN(i), AA(i)]);
end

figure;
surf(NN, log10(AA), Z);
xlabel('n');
ylabel('log_{10}(A)');
zlabel('Function Value');
title('Plot of the Function');