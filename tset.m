%% Assignment
% Differentiation

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;
close all;


f = @(x) cos(x); 
f_actual = @(x) -sin(x);

x = 2;


function [val] = forward_method(f, x, h)
    val = (f(x+h) - f(x)) / h;
end

function [val] = backward_method(f, x, h)
    val = (f(x) - f(x-h)) / h;
end

function [val] = central_method(f, x, h)
    val = (f(x+h) - f(x-h)) / (2*h);
end

function [val] = central_method_h4(f, x, h)
    val = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
end


steps = 50;
hs = logspace(-2, -10, steps);



fx            = f(x) * ones(1, steps);          
fx_h          = f(x) * ones(1, steps);  


fx           = zeros(1, steps);
fx_derivative = zeros(1, steps);

for i = 1:steps
    fx(i)   = f(x);
    fx_h(i) = f(x+hs(i));
end


f_d = zeros(1, steps);
err = zeros(1, steps);

for i = 1:steps
    f_d(i) = (fx_h(i) - fx(i))/hs(i);
    err(i) = abs(f_d(i) - f_actual(x));
end

loglog(hs, err, 'DisplayName', 'Forward error')
hold on




