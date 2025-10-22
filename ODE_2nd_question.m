clc;
clear;
close all;

% Parameters
alpha = 1;    
h = 0.1;      
k = 0.001;    

x_start = 0;
x_end = 1;
t_start = 0;
t_end = 0.1;

% Discretize domain
xs = x_start:h:x_end;
ts = t_start:k:t_end;

nx = length(xs);
nt = length(ts);


T = zeros(nx, nt);

for i = 1:nx
    T(i,1) = sin(pi * xs(i)) + sin(2*pi * xs(i));
end

T(1,:) = 0;
T(end,:) = 0;

for j = 1:nt-1
    for i = 2:nx-1
        T(i,j+1) = T(i,j) + (alpha * k / (h^2)) * (T(i+1,j) - 2*T(i,j) + T(i-1,j));
    end
end

[X, Y] = meshgrid(xs, ts);

figure;
contourf(X, Y, T');

figure;
surf(T)
colorbar;
