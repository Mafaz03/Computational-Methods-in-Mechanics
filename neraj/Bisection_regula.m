%% Assignemnt - 4
% Bisection and Regula Falsi

% Name        : Niraj Kumar Singh
% Roll Number : AM25M807
% Department  : Applied Mechanics

%% Part 1 (Preprocessing)
clc; clear; close all;
f  = @(x) sin(10*x) + cos(3*x);
tol = 1e-4;

%% Part 2 (Processing)

% Bisection for root near 3.74575
a = 3;  b = 6;  r1 = 3.74575;
if f(a)*f(b) > 0, error('Bad bracket'); end

k = 0;  it1 = []; err1 = [];
c = (a+b)/2;
while abs(f(c)) > tol
    if f(a)*f(c) < 0, b = c; else, a = c; end
    c = (a+b)/2;  k = k + 1;
    it1(k)  = k;
    err1(k) = abs(c - r1);
    fprintf('Bisection %2d  c=%.6f\n', k, c);
end

semilogy(it1, err1, 'r-o','LineWidth',1.2); hold on
xlabel('Iteration'); ylabel('|c - root|')
title('True Error vs Iteration')

R = arrayfun(@(j) log(err1(j)/err1(j-1))/ ...
                  log(err1(j-1)/err1(j-2)), 3:numel(err1));
fprintf('\nBisection order ≈ %.4f\n\n', mean(R));


% Regula-Falsi for root near 5.67903
a = 3;  b = 6;  r2 = 5.67903;
if f(a)*f(b) > 0, error('Bad bracket'); end

k = 0;  it2 = []; err2 = [];
c = b - f(b)*(b - a)/(f(b) - f(a));
while abs(f(c)) > tol
    if f(a)*f(c) < 0, b = c; else, a = c; end
    c = b - f(b)*(b - a)/(f(b) - f(a));
    k = k + 1;
    it2(k)  = k;
    err2(k) = abs(c - r2);
    fprintf('Regula-Falsi %2d  c=%.6f\n', k, c);
end

%% Part 3 (Post Processing / Plotting)

semilogy(it2, err2, 'b-s','LineWidth',1.2)
legend('Bisection','Regula-Falsi','Location','northeast')

R = arrayfun(@(j) log(err2(j)/err2(j-1))/ ...
                  log(err2(j-1)/err2(j-2)), 3:numel(err2));
fprintf('\nRegula-Falsi order ≈ %.4f\n', mean(R));