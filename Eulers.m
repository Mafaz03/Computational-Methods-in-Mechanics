clc;
clear;
close all;

%% Eulers
f =@(y, x) (-1.2 * y) + (7 * exp(-0.3*x));
h = 0.1;
a = 0;
b = 3;
y = 3; % Initial Value

M = (b-a)/h;
% x = linspace(a,b,M);
x = a:h:b;

ys = [];
for i = 1: M
    y = y + (h * f(y, x(i)));
    ys(end+1) = y;
end
disp(table(ys', 'VariableNames', {'Euler'}))

figure()
plot(x(1:end-1), ys)
title("Tracing x | Euler")
ylabel('x')
xlabel('time')

%% RK 2
f   =@(y, x) y + x^3;
f_a =@(y, x) (7 * exp(x)) - (x^3) - (3 * x^2) - (6 * x^2) - 6;
h = 0.1;
a = 0;
b = 1.5;
y = 1; % Initial Value
y_a = 1;
M = (b-a)/h;
% x = linspace(a,b,M);
x = a:h:b;
fprintf("\n\nRK 2\n")
fprintf("y        | Error\n")
ys = [];
for i = 1: M
    k1 = h * f(y, x(i));
    k2 = h * f(y + k1/2, x(i) + (h/2));
    y  = y + k2;
    ys(end+1) = y;

    y_a = f_a(y_a, i*h);
    fprintf("%f | %f\n", y, abs(y - y_a))
end
% disp(table(ys', 'VariableNames', {'RK 2'}))

figure()
plot(x(1:end-1), ys)
title("Tracing x | RK 2")
ylabel('x')
xlabel('time')

%% RK 4
f   =@(y, x) y + x^3;
f_a =@(y, x) (7 * exp(x)) - (x^3) - (3 * x^2) - (6 * x^2) - 6;
h = 0.1;
a = 0;
b = 1.5;
y = 1; % Initial Value
y_a = 1;
M = (b-a)/h;
% x = linspace(a,b,M);
x = a:h:b;

fprintf("\n\nRK 4\n")
fprintf("y        | Error\n")
ys = [];
for i = 1: M
    k1 = h * f(y, x(i));
    k2 = h * f(y + k1/2, x(i) + (h/2));
    k3 = h * f(y + k2/2, x(i) + (h/2));
    k4 = h * f(y + k3, x(i) + h);
    y = y + (k1/6) + (k2/3) + (k3/3) + (k4/6);
    ys(end+1) = y;

    y_a = f_a(y_a, i*h);
    fprintf("%f | %f\n", y, abs(y - y_a))
end

figure()
plot(x(1:end-1), ys)
title("Tracing x | RK 4")
ylabel('x')
xlabel('time')