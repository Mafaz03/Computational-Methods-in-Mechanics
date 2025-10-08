clc;
clear;
close all;

global sigma rho beta

% Parameters
sigma = 10; 
rho   = 28;   
beta  = 8/3;  

% Lorenz equations
fx_t = @(x, y) sigma * (y - x);
fy_t = @(x, y, z) x * (rho - z) - y;
fz_t = @(x, y, z) (x * y) - (beta * z);

% Initial conditions
x_ini = 1;
y_ini = 1;
z_ini = 1;

% Time setup
t_final = 50;
h = 0.01;
t = 0:h:t_final;
N = length(t);

% Allocate
x = zeros(1, N);
y = zeros(1, N);
z = zeros(1, N);

x(1) = x_ini;
y(1) = y_ini;
z(1) = z_ini;

%% Euler
for i = 1:N-1
    x(i+1) = x(i) + h * fx_t(x(i), y(i));
    y(i+1) = y(i) + h * fy_t(x(i), y(i), z(i));
    z(i+1) = z(i) + h * fz_t(x(i), y(i), z(i));
end

figure;
plot3(x, y, z);
title('Lorenz System - Euler Method');
xlabel('x'); ylabel('y'); zlabel('z');
grid on;

hold on

%% RK4
x_rk(1) = x_ini;
y_rk(1) = y_ini;
z_rk(1) = z_ini;

for i = 1:N-1
    % k1
    k1x = fx_t(x(i), y(i));
    k1y = fy_t(x(i), y(i), z(i));
    k1z = fz_t(x(i), y(i), z(i));

    % k2
    k2x = fx_t(x(i) + 0.5*h*k1x, y(i) + 0.5*h*k1y);
    k2y = fy_t(x(i) + 0.5*h*k1x, y(i) + 0.5*h*k1y, z(i) + 0.5*h*k1z);
    k2z = fz_t(x(i) + 0.5*h*k1x, y(i) + 0.5*h*k1y, z(i) + 0.5*h*k1z);

    % k3
    k3x = fx_t(x(i) + 0.5*h*k2x, y(i) + 0.5*h*k2y);
    k3y = fy_t(x(i) + 0.5*h*k2x, y(i) + 0.5*h*k2y, z(i) + 0.5*h*k2z);
    k3z = fz_t(x(i) + 0.5*h*k2x, y(i) + 0.5*h*k2y, z(i) + 0.5*h*k2z);

    % k4
    k4x = fx_t(x(i) + h*k3x, y(i) + h*k3y);
    k4y = fy_t(x(i) + h*k3x, y(i) + h*k3y, z(i) + h*k3z);
    k4z = fz_t(x(i) + h*k3x, y(i) + h*k3y, z(i) + h*k3z);

    % Update
    x_rk(i+1) = x_rk(i) + (h/6)*(k1x + 2*k2x + 2*k3x + k4x);
    y_rk(i+1) = y_rk(i) + (h/6)*(k1y + 2*k2y + 2*k3y + k4y);
    z_rk(i+1) = z_rk(i) + (h/6)*(k1z + 2*k2z + 2*k3z + k4z);
end

plot3(x_rk, y_rk, z_rk);
title('Lorenz System - RK4 Method');
xlabel('x'); ylabel('y'); zlabel('z');
grid on;

legend("Euler", "RK4")


errors = zeros(1, N);
for i = 1:N
    errors(i) = ((x(i) - x_rk(i))^2 + (y(i) - y_rk(i))^2 + (z(i) - z_rk(i))^2)^0.5;
end

figure
plot(t, errors)
xlabel("Time")
ylabel("Errors")
title("Error RK4 & Euler")