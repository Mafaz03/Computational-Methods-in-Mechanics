%% Assignemnt - 4
% Newton and Secant

% Name        : Niraj Kumar Singh
% Roll Number : AM25M807
% Department  : Applied Mechanics

%% Part 1 (Preprocessing)

clc; clear; close all;

% --------- Secant Method ----------
x0 = 10; x1 = 9; tol = 1e-6;
f = @(x) (x-1).*(exp(x-1)-1);

rel_err = inf; k = 0;
it_s = []; err_s = [];

%% Part 2 (Processing)
old_err = 0;

fprintf('++++++ SECANT METHOD ++++++\n\n')
while rel_err > tol
    x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0)) ;
    rel_err = abs((x2 - x1)/x2);
    x0 = x1; x1 = x2;
    k = k + 1;
    it_s(k)  = k;
    err_s(k) = rel_err;
    fprintf('Loop %d  |  c = %.6f\n', k, x2);
end

semilogy(it_s, err_s, 'r-s','LineWidth',1.2); hold on
xlabel('Iteration'); ylabel('|c - root|')
title('True Error vs Iteration')

R = arrayfun(@(j) log(err_s(j)/err_s(j-1))/ ...
                  log(err_s(j-1)/err_s(j-2)), 3:numel(err_s));
fprintf('\nEstimated order (Secant): %.4f\n\n', mean(R));


% --------- Newton Method ----------
x0 = 10; tol = 1e-6;
syms x
df = matlabFunction(diff((x-1)*(exp(x-1)-1)));

rel_err = inf; k = 0;
it_n = []; err_n = [];

fprintf('++++++ NEWTON METHOD ++++++\n\n')
while rel_err > tol
    x1 = x0 - f(x0)/df(x0);
    rel_err = abs((x1 - x0)/x1);
    x0 = x1;
    k = k + 1;
    it_n(k)  = k;
    err_n(k) = rel_err;
    fprintf('Loop %d  |  c = %.6f\n', k, x1);
end

%% Part 3 (Post Processing / Plotting)

semilogy(it_n, err_n, 'y-o','LineWidth',1.2)
legend('Secant','Newton','Location','northeast')

R = arrayfun(@(j) log(err_n(j)/err_n(j-1))/ ...
                  log(err_n(j-1)/err_n(j-2)), 3:numel(err_n));
fprintf('\nEstimated order (Newton): %.4f\n', mean(R));