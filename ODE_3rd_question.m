clc;
clear;
close all;

h = 0.1;
k = 0.1;

x_start = 0;
x_end = 1;
t_start = 0;
t_end = 1;

xs = x_start:h:x_end;
ts = t_start:k:t_end;

top_bc = 0;
bottom_bc = 0;

matrix = zeros(length(xs), length(ts));
matrix(1, :) = top_bc;
matrix(end, :) = bottom_bc;


for i = 1:length(xs)
    matrix(i, 1) = sin(pi * xs(i));  
end

r = (k / h);
for i = 2:length(xs)-1
    matrix(i, 2) = matrix(i, 1) + (r^2)/2 * (matrix(i-1, 1) - 2 * matrix(i, 1) + matrix(i+1, 1));
end

for j = 2:length(ts)-1
    for i = 2:length(xs)-1

        matrix(i, j+1) = 2*matrix(i, j) - matrix(i, j-1) + r^2*(matrix(i+1, j) - 2*matrix(i, j) + matrix(i-1, j));
    end
end


[X, Y] = meshgrid(xs, ts);
figure;
contourf(X, Y, matrix');
title('Contour lines');
colorbar;

figure;
surf(matrix)
colorbar;
