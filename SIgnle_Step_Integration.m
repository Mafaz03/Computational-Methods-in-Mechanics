clc;
clear;
close all;

answer = 0.74306944;

% f =@(x) exp(-x^2);
f =@(x) 1 / ((1+x^3)^0.5);

a = 1;
b = 3;

%% Single Step

% Trapezoid
h = (b-a);
trap_integral = h/2 * (f(a) + f(b));
fprintf("Trapezoid Single step: %f\n\n", trap_integral)

% Simpsion's 1/3rd
h = (b-a)/2;
simp_1_3_integral = (b-a)/6 * (f(a) + 4*f(a+h) + f(b));
fprintf("Simpsion's 1/3rd Single step: %f\n", simp_1_3_integral)

% Simpsion's 3/8rd
h = (b-a)/3;
simp_3_8_integral = (3*h/8) * (f(a) + 3*f(a+h) + 3*f(a+2*h) + f(b));
fprintf("Simpsion's 3/8th Single step: %f\n", simp_3_8_integral)


%% Multiple Steps

% Trapezoid
trap_mul_hs = [];
trap_mul_errors = [];
M = linspace(100, 200, 10);

% M = 100;
for m = 1:length(M)
    h = (b-a)/m;
    trap_mul_hs(end+1) = h;
    k_val = linspace(a, b, m);
    
    sum = 0;
    for k = 1:m-1
        sum = sum + (   f(k_val(k))  +  f(k_val(k+1))   );
    end
    trap_mul = (h/2) * sum;
    trap_mul_errors(end+1) = abs(answer - trap_mul);
end
figure()
plot(trap_mul_hs, trap_mul_errors)
xlabel('h');
ylabel('abs error');
title("Trapezoid Multi-Step optimal h value")
fprintf("\n\nTrapezoid Multi step: %f\n",trap_mul)


% Simpsion's 1/3rd

simp_1_3_mul_hs = [];
simp_1_3_mul_errors = [];
M = linspace(100, 200, 10);
for m = 1:length(M)
    n = 2 * m; % n = 2M
    k_val = linspace(a, b, n/2);
    h = (b-a)/m;
    simp_1_3_mul_hs(end+1) = h;

    sum = 0;
    
    for k = 1:n/2
        if mod(k, 2) == 0
            sum = sum + (2 * f(k_val(k)));
        else
            sum = sum + (4 * f(k_val(k)));
        end
    end
    
    simp_1_3_integral = (h/3) * (f(a) + sum + f(b));
    simp_1_3_mul_errors(end+1) = abs(simp_1_3_integral-answer);
end
fprintf("Simpsion's 1/3rd Multi step: %f\n\n\n", simp_1_3_integral)

figure()
plot(simp_1_3_mul_hs, simp_1_3_mul_errors)
xlabel('h');
ylabel('abs error');
title("Simpsion's 1/3rd Multi-Step optimal h value")

%% Gauss Legendre

f_gl =@(x) 1 / sqrt(1+((x+2)^3));

% 2 point
x_2_1 = -1/sqrt(3);
x_2_2 = 1/sqrt(3);

w_2_1 = 1;
w_2_2 = 1;

gl_2 = (w_2_1 * f_gl(x_2_1)) + (w_2_2 * f_gl(x_2_2));
fprintf("Gauss Legendre 2 point: %f\n", gl_2)

% 3 point
x_3_1 = -sqrt(3/5);
x_3_2 = 0;
x_3_3 = sqrt(3/5);

w_3_1 = 5/9;
w_3_2 = 8/9;
w_3_3 = 5/9;

gl_3 = (w_3_1 * f_gl(x_3_1)) + (w_3_2 * f_gl(x_3_2)) + (w_3_3 * f_gl(x_3_3));
fprintf("Gauss Legendre 3 point: %f\n\n", gl_3)

% Adaptive Quadrature using Simpson's 1/3 Multi step method
tol = 1e-5;

function [simp_1_3_integral, h] = Simposon1_3 (f, a, b)
    M = 200;
    n = 2 * M; % n = 2M
    k_val = linspace(a, b, n/2);
    h = (b-a)/M;
    sum = 0;

    for k = 1:n/2
        if mod(k, 2) == 0
            sum = sum + (2 * f(k_val(k)));
        else
            sum = sum + (4 * f(k_val(k)));
        end
    end

    simp_1_3_integral = (h/3) * (f(a) + sum + f(b));
end

global operations hs
operations = 0;
hs = [];

function I = AdaptiveSimpson(f, a, b, tol)
    global operations hs

    operations = operations + 1;
    c = (a+b)/2;
    [S1, ~] = Simposon1_3(f, a, c);
    [S2, ~] = Simposon1_3(f, c, b);

    [I1, h] = Simposon1_3(f, a, b);
    I2 = S1 + S2;

    hs(end+1) = h;

    if abs(I1 - I2) < 15 * tol
        I = I2 + (I2 - I1)/15;

    else
        I = AdaptiveSimpson(f, a, c, tol/2) + ...
            AdaptiveSimpson(f, c, b, tol/2);
    end
end

A_Simposn = AdaptiveSimpson(f, a, b, tol);
fprintf("Adaptive Quadrature using Simpson: %f, number of operations: %d\n\n", A_Simposn, operations)
fprintf("h that were changed in adaptive simpson:\n")
disp(table(hs', 'VariableNames', {'h'}))