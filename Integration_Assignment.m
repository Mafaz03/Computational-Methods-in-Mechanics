clc;
clear;
close all;

f = @(x) 1 ./ sqrt(1 + x.^3);

a = 1;
b = 3;
M = 100;

I_exact = integral(f, a, b);

%% Trapezoidal Rule
h = (b - a) / M;
x = linspace(a, b, M+1); 
fx = f(x);
trap_mul = (h/2) * (fx(1) + 2*sum(fx(2:end-1)) + fx(end));
err_trap = abs(I_exact - trap_mul);

fprintf("Trapezoid Multi step: %f\n",trap_mul)
fprintf("Absolute Error (Trapezoidal): %e\n\n", err_trap);

%% Simpson's 1/3 Rule
n = 2*M; 
h = (b-a)/n;
x = linspace(a, b, n+1);
fx = f(x);

% Apply Simpson's coefficients: 1, 4, 2, 4, 2,...,4, 1
coeff = ones(1, n+1); 
coeff(2:2:n) = 4;     
coeff(3:2:n-1) = 2;   

simp_1_3_integral = (h/3) * sum(coeff .* fx);
err_simp = abs(I_exact - simp_1_3_integral);

fprintf("Simpson's 1/3 Multi step: %f\n", simp_1_3_integral)
fprintf("Absolute Error (Simpson's 1/3): %e\n", err_simp);
