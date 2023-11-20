clear, clc;

data = load("Sharad.mat");
A = 10^-23;
a = data.a; dhdx = data.dhdx; g = data.g; rho = data.rho; x = data.x; H_obs = data.H_obs;

H = @(n, A) (-(2+n)/(2*A).*a.*(rho*g).^(-n).*abs(dhdx).^(1-n).*dhdx.^-1).^(1./(n+2));

sum1 = @(n, A) sum((H(n,A)- H_obs).^2);

fun = @(args) sqrt(sum1(args(1),args(2)));
initialdata = [2.89, 2.09e-26];



options = optimoptions('fminunc');
options.TolFun = 1e-35;
options.OptimalityTolerance = 1e-35;    % This completely sucks
options.StepTolerance = 1e-35;
options.MaxFunctionEvaluations = 5e8;
options.MaxIterations = 1e8;
[minParams, fval] = fminunc(fun, initialdata, options);




%%


n_values = linspace(2.5, 4, 1000);
A_values = linspace(0.01e-25, 0.00001e-25, 1000);

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

% figure;
% surf(NN, AA, Z);
% xlabel('n');
% ylabel('log_{10}(A)');
% zlabel('Function Value');
% title('Plot of the Function');

