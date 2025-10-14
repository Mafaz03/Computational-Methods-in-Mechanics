clc
clear
close all

nodes = 50;

matrix      = zeros(nodes,nodes);

% Boundary Condition
matrix(1, :)     = 100; % Bottom
matrix(nodes, :) = 10;  % Top

matrix(:, 1)     = 0;   % Left
matrix(:, nodes) = 0;   % Right


tol = 1e-5;
abs_err = tol+1;
loops = 1;

while abs_err > tol
    matrix_old = matrix;
    for i = 2:nodes-1
        for j = 2:nodes-1
           
            matrix(i,j) = (matrix_old(i+1, j) + matrix_old(i-1, j) + matrix_old(i, j+1) + matrix_old(i, j-1)) / 4;
        end
    end

    abs_err = max(max(abs(matrix - matrix_old)));
    loops = loops+1;
    matrix_old = matrix;
    loops
end

[X, Y] = meshgrid(1:nodes, 1:nodes);

% Contour lines
figure;
contour(X, Y, matrix);
title('Contour lines');

% Filled contour
figure;
contourf(X, Y, matrix);
title('Filled contour');
colorbar;
