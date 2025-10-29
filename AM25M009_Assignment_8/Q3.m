%% Assignemnt
% PDE's

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;
close all;

%% Part 1 (Preprocessing)
h = 0.01;
k = 0.01;

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

%% Part 2 (Processing / Using the Algorithms)
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
title(sprintf('Contour lines for r = %f', r));
colorbar;

figure;
surf(matrix)
title(sprintf('Surf for r = %f', r));
colorbar;

h = 0.01;
k = 0.01;

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

    plot(xs, matrix(:, j+1));
    title(sprintf('Time = %f', ts(j+1)));
    grid on;
    drawnow;
    pause(0.01)
end

%% Part 3 (Post Processing / Plotting)
% Final plots
[X, Y] = meshgrid(xs, ts);
figure;
contourf(X, Y, matrix');
title('Contour lines');
colorbar;

figure;
surf(matrix)
colorbar;