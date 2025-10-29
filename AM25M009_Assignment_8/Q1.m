%% Assignemnt
% PDE's

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;
close all;


%% Part 1 (Preprocessing)
% Parameters
width = 1.0;      
height = 1.5;     
k = 0.4;          
Q_dot = 100;     
nx = 5;          
ny = 7;    

max_iterations = 100;

% Grid and Coefficients
dx = width / (nx - 1);
dy = height / (ny - 1);
alpha = 1 / (dx^2);
beta = 1 / (dy^2);
source = -Q_dot / k;

omega = 0.1; 

T = zeros(ny,nx);
tolerance = 1e-2;

tic;
%% Part 2 (Processing / Using the Algorithms)
omegas = [0.25, 0.5, 1];
for p = 1:length(omegas)
    omega = omegas(p)
    for iter = 1:max_iterations
        T_old = T;
    
        % Interior points
        for j = 2:ny-1
            for i = 2:nx-1
                T_new = (alpha * (T(j, i+1) + T(j, i-1)) + ...
                         beta  * (T(j+1, i) + T(j-1, i)) - ...
                         source) / (2*alpha + 2*beta);
    
                % relaxation
                T(j,i) = (1 - omega) * T(j,i) + omega * T_new;
            end
        end
    
        % Neumann boundaries
        for j = 2:ny-1
            % Left
            T_new_left = (2 * alpha * T(j, 2) + ...
                          beta * (T(j+1, 1) + T(j-1, 1)) - ...
                          source) / (2*alpha + 2*beta);
            T(j,1) = (1 - omega) * T(j,1) + omega * T_new_left;
    
            % Right
            T_new_right = (2 * alpha * T(j, nx-1) + ...
                           beta * (T(j+1, nx) + T(j-1, nx)) - ...
                           source) / (2*alpha + 2*beta);
            T(j,nx) = (1 - omega) * T(j,nx) + omega * T_new_right;
        end
    
    
        max_change = max(max(abs(T - T_old)));
        if max_change < tolerance
            fprintf('Converged after %d iterations\n', iter);
            break;
        end
    end
    
    elapsed_time = toc;
    fprintf('Elapsed time: %f seconds\n for 5 points', elapsed_time);
    
    
    disp('Temperature Distribution:');
    disp(T);
    
    figure;
    surf(T)
    title(sprintf("r = %f", omega))
    
    x_coords = linspace(0, width, nx);
    y_coords = linspace(0, height, ny);
    
    [X, Y] = meshgrid(x_coords, y_coords);
    figure;
    contourf(X, Y, T);
    title(sprintf("r = %f", omega))
    hold on;
    % [C, h] = contour(X, Y, T);
    % colorbar;
    % xlabel('x');
    % ylabel('y');
end


Width  = 1.0;
Height = 1.5;
alpha  = 1;
beta   = 1;
source = -source;
tolerance = 1e-6;
max_iterations = 10000;
omega = 0.1; % relaxation factor

T = zeros(ny, nx);

tic;
for iter = 1:max_iterations
    T_old = T;
    max_change = 0;

    % Interior points
    for j = 2:ny-1
        for i = 2:nx-1
            % 9-point 
            T_new = ( ...
                4*(T(j, i+1) + T(j, i-1) + T(j+1, i) + T(j-1, i)) + ...
                (T(j+1, i+1) + T(j+1, i-1) + T(j-1, i+1) + T(j-1, i-1)) + ...
                6*source*dx^2 ) / 20;

            % relaxation
            T(j,i) = T(j,i) + omega * (T_new - T(j,i));

            % Track error
            err = abs(T(j,i) - T_old(j,i));
            if err > max_change
                max_change = err;
            end
        end
    end

    % Neumann boundaries
    for j = 2:ny-1
        % Left
        T_new_left = T(j,2);
        T(j,1) = (1 - omega) * T(j,1) + omega * T_new_left;

        % Right
        T_new_right = T(j,nx-1);
        T(j,nx) = (1 - omega) * T(j,nx) + omega * T_new_right;
    end

    % Convergence check
    if max_change < tolerance
        break;
    end
end

elapsed_time = toc;
fprintf('Elapsed time: %f seconds\n for 9 points', elapsed_time);
T

%% Part 3 (Post Processing / Plotting)

figure;
surf(T);
title('Temperature Distribution Crank Nicholson');
xlabel('x'); ylabel('y'); zlabel('T');

x_coords = linspace(0, Width, nx);
y_coords = linspace(0, Height, ny);
[X, Y] = meshgrid(x_coords, y_coords);

figure;
contourf(X, Y, T);
title('Temperature Distribution Crank Nicholson');
hold on;
colorbar;








