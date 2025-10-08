clc;
clear;
close all;

global sigma rho beta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 10; 
rho   = 28;   
beta  = 8/3;  

fx_t =@(x, y)    sigma * (y-x);
fy_t =@(x, y, z) (x * (rho - z)) - y;
fz_t =@(x, y, z) (x*y) - (beta*z);

x_ini = 1;
y_ini = 1;
z_ini = 1;

t_final = 50;
h = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:h:t_final;

x(1) = x_ini;
y(1) = y_ini;
z(1) = z_ini;

%% Eulers
y_exact = zeros(size(t));
for i = 1:length(t)
    x(i+1) = x(i) + h * fx_t(x(i), y(i));
    y(i+1) = y(i) + h * fy_t(x(i), y(i), z(i));
    z(i+1) = z(i) + h * fz_t(x(i), y(i), z(i));
end

plot3(x,y,z)

%% RK 4
x(1) = x_ini;
y(1) = y_ini;
z(1) = z_ini;

% for i = 1: length(t)
%     k1 = h * f(y, x(i));
%     k2 = h * f(y + k1/2, x(i) + (h/2));
%     k3 = h * f(y + k2/2, x(i) + (h/2));
%     k4 = h * f(y + k3, x(i) + h);
%     y = y + (k1/6) + (k2/3) + (k3/3) + (k4/6);
% 
%     ys(end+1) = y;
% 
% end



