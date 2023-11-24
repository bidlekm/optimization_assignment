clear, clc;

data = load("Sharad.mat");
A = 10^-23;
a = data.a; dhdx = data.dhdx; g = data.g; rho = data.rho; x = data.x; H_obs = data.H_obs;

H = @(n, A) (-(2+n)/(2*10.^A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^-1).^(1./(n+2));

sum1 = @(n, A) sum((H(n,A)- H_obs).^2);

fun = @(args) sqrt(sum1(args(1),args(2)));
initialdata = [2.89, -25];



options = optimoptions('fminunc');
% options.TolFun = 1e-35;
% options.OptimalityTolerance = 1e-35;    % This completely sucks
% options.StepTolerance = 1e-35;
% options.MaxFunctionEvaluations = 5e8;
% options.MaxIterations = 1e8;
[minParams, fval] = fminunc(fun, initialdata, options);




%%


n_values = linspace(2, 10, 100);
A_values = linspace(-100, 0, 100);

[NN, AA] = meshgrid(n_values, A_values);
Z = zeros(size(NN));

for i = 1:numel(NN)
    Z(i) = fun([NN(i), AA(i)]);
end

zmin = min(Z(:));
[row, col] = find(Z == zmin);
nmin = n_values(col)
Amin = A_values(row)
fun([nmin, Amin])

% Real values
% Amin = 1.5051e-26, nmin = 2.90394

figure;
surf(NN, AA, log10(Z));
xlabel('n');
ylabel('log_{10}(A)');
zlabel('Function Value');
title('Plot of the Function');

