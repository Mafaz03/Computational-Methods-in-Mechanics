%% Assignemnt
% PDE's

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

%% Part 1 (Preprocessing)
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

% Initialize
T = zeros(nx, nt);       % FTCS
T_CN = zeros(nx, nt);    % Crank–Nicolson

% Initial condition
for i = 1:nx
    T(i,1) = sin(pi * xs(i)) + sin(2*pi * xs(i));
    T_CN(i,1) = T(i,1);
end

T(1,:) = 0; 
T(end,:) = 0;
T_CN(1,:) = 0; 
T_CN(end,:) = 0;


%% Part 2 (Processing / Using the Algorithms)
% FTCS
for j = 1:nt-1
    for i = 2:nx-1
        T(i,j+1) = T(i,j) + (alpha * k / h^2) * (T(i+1,j) - 2*T(i,j) + T(i-1,j));
    end
end

%% Crank–Nicolson 
r = alpha * k / (2 * h^2);
A = diag((2 + 2*r) * ones(nx, 1)) + diag(-r * ones(nx-1, 1), 1) + diag(-r * ones(nx-1, 1), -1);

% Boundary conditions
A(1,:) = 0; 
A(1,1) = 1;
A(nx,:) = 0; 
A(nx,nx) = 1;

for j = 1:nt-1
    Tj = T_CN(:, j);
    d = zeros(nx, 1);
    d(1) = 0;
    d(nx) = 0;
    for i = 2:nx-1
        d(i) = r*Tj(i-1) + (2 - 2*r)*Tj(i) + r*Tj(i+1);
    end
    T_CN(:, j+1) = A \ d;
end


% r = alpha * k / (2 * h^2);
% A = diag((2 + 2*r) * ones(nx, 1)) + diag(-r * ones(nx-1, 1), 1) + diag(-r * ones(nx-1, 1), -1);
% 
% Boundary conditions
% A(1,:) = 0; 
% A(1,1) = 1;
% A(nx,:) = 0; 
% A(nx,nx) = 1;
% 
% tolerance = 1e-10;      
% max_iter = 1000;        
% 
% T_history = [];
% 
% for j = 1:nt-1
%     Tj = T_CN(:, j);
%     d = zeros(nx, 1);
%     d(1) = 0; d(nx) = 0;
%     for i = 2:nx-1
%         d(i) = r*Tj(i-1) + (2 - 2*r)*Tj(i) + r*Tj(i+1);
%     end
% 
%     T_new = T_CN(:, j);  % initial guess
%     T_history = zeros(max_iter, nx);  
% 
%     for iter = 1:max_iter
%         T_old = T_new;
%         for i = 2:nx-1
%             T_new(i) = (r*T_new(i-1) + (2 - 2*r)*Tj(i) + r*T_new(i+1)) / (2 + 2*r);
%         end
%         T_history(iter, :) = T_new;  
%         if max(abs(T_new - T_old)) < tolerance
%             T_history = T_history(1:iter, :);  
%             break;
%         end
%     end
% 
%     T_CN(:, j+1) = T_new;
% end

%% Part 3 (Post Processing / Plotting)
[X, Y] = meshgrid(xs, ts);

figure;
contourf(X, Y, T');
title('FTCS (Explicit)');
xlabel('x'); ylabel('t'); colorbar

figure;
surf(T);
title('FTCS (Explicit)');
xlabel('x'); ylabel('t'); colorbar

figure;
contourf(X, Y, T_CN');
title('Crank–Nicolson (Implicit)');
xlabel('x'); ylabel('t'); colorbar


figure;
surf(T_CN);
title('Crank–Nicolson (Implicit)');
xlabel('x'); ylabel('t'); colorbar
