%% Assignment
% Differentation

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
    val = (f(x+h) - f(x))/h;
end

function [val] = backward_method(f, x, h)
    val = (f(x) - f(x-h))/ h ;
end

function [val] = central_method(f, x, h)
    val = ( f(x+h) - f(x-h) ) / (2*h) ;
end

function [val] = central_method_h4(f, x, h)
    val = ( -f(x+2*h) + (8 * f(x+h)) - (8 * f(x-h)) + f(x-2*h)) / (12*h) ;
end


steps = 50;
hs = logspace(-2, -10, steps);

fx = zeros(1, steps);
for i = 1:steps
    fx(i) = f(hs(i));
end


error_forward = zeros(1, steps);
error_backward = zeros(1, steps);
error_central = zeros(1, steps);
error_central_h4 = zeros(1, steps);

for i = 1:steps
    forwardVal  = forward_method(f, x, hs(i));
    backwardVal = backward_method(f, x, hs(i));
    centralVal = central_method(f, x, hs(i));
    centralVal_h4 = central_method_h4(f, x, hs(i));

    error_forward(i)  = abs(forwardVal - f_actual(x));
    error_backward(i) = abs(backwardVal - f_actual(x));
    error_central(i)  = abs(centralVal - f_actual(x));
    error_central_h4(i)  = abs(centralVal_h4 - f_actual(x));
end



loglog(hs, error_forward, 'DisplayName', 'Forward error')
hold on
loglog(hs, error_backward, 'DisplayName', 'Backward error')
hold on
loglog(hs, error_central, 'DisplayName', 'Central error')
hold on
loglog(hs, error_central_h4, 'DisplayName', 'Central error O(h^4)')

[~, min_h_a] = min(error_forward);
[~, min_h_b] = min(error_backward);
[~, min_h_c] = min(error_central);
[~, min_h_d] = min(error_central_h4);

fprintf("O(h) Forward hopt error: %d min error: %d\n", hs(min_h_a), min(error_forward));
fprintf("O(h) Backward hopt error: %d min error: %d\n", hs(min_h_b), min(error_backward))
fprintf("O(h^2) Central hopt error: %d min error: %d\n", hs(min_h_c), min(error_central))
fprintf("O(h^4) Central hopt error: %d min error: %d\n", hs(min_h_d), min(error_central_h4))



legend show;
% grid on;
xlabel('h');
ylabel('abs error');

