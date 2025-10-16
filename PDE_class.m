% clc
% clear
close all

nodes = 10;

matrix      = zeros(nodes,nodes);

% Boundary Condition
matrix(1, :)     = 100; % Bottom
matrix(nodes, :) = 10;  % Top

matrix(:, 1)     = 0;   % Left
matrix(:, nodes) = 0;   % Right

coses = cos(pi/(nodes-1)) + cos(pi/(nodes-1));
w = 4 / (2 + sqrt(4 - coses^2)); % change this later on
% w = 1;

tol = 1e-5;
abs_err = tol+1; % To allow it once lol
loops = 1;

while abs_err > tol
    matrix_old = matrix;
    for i = 2:nodes-1
        for j = 2:nodes-1
           
            % No relaxation
            % matrix(i,j) = (matrix(i+1, j) + matrix(i-1, j) + matrix(i, j+1) + matrix(i, j-1)) / 4;

            % relaxation
            s =  (matrix_old(i+1, j) + matrix(i-1, j) + matrix(i, j+1) + matrix(i, j-1) - 4 * matrix_old(i, j)) / 4;
            matrix(i,j) = matrix_old(i, j) + (w * s);

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

% contour
figure;
contourf(X, Y, matrix);
title('Filled contour');
colorbar;
