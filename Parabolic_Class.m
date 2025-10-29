clc;
clear;
close all;


m = 10;
dx = 2;
x  = linspace(0, dx, m);
t = 10;
c = 0.12;

r = c / dx^2;

u = zeros(m, t);

u(1, :) = 0;      % BC
u(end, :) = 10;   % BC

for i = 1:m
    if i <= m/2
        u(i, 1) = 100 * x(i);
    else
        u(i, 1) = 200 - (100 * x(i));
    end
end

for j = 1:t-1
    for i = 2:m-1
        u(i, j+1) = u(i, j) + r * (u(i+1, j) - 2*u(i, j) + u(i-1, j));
    end
end

surf(u)
