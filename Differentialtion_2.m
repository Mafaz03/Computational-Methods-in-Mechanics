clc;
clear;
close all;

syms x
expression = exp(x) * sin(x);

% function and derivative as numeric function handles
f        = matlabFunction(expression);
f_actual = matlabFunction(diff(expression, x));



forward_method  = @(f, x, h) (f(x+h) - f(x)) / h;
backward_method = @(f, x, h) (f(x) - f(x-h)) / h;
central_method  = @(f, x, h) (f(x+h) - f(x-h)) / (2*h);

methods = {forward_method, backward_method, central_method};
methods_name = ["Forward Method", "Backward Method", "Central Method"];

function [] = Derivative(method_name, func, f, f_actual, val, h)

    h_values = [];
    relative_errors = [];

    fprintf("\n\n========%s========\n\n", method_name)
    tolerance = 1e-6;
    loop = 0;
    relative_error = 0;
    while (relative_error > tolerance) || loop == 0
        pred_val  = func(f, val, h);
        actual_val   = f_actual(val);
        
        relative_error = abs((pred_val - actual_val) / actual_val);

        h_values(end+1) = h;
        relative_errors(end+1) = relative_error;

        h    = h/1.5;
        loop = loop + 1;
        
        fprintf("%d | %.8f | %.8f | error = %.2e\n", loop, pred_val, actual_val, relative_error);
    end
   
    plot(relative_errors,h_values)
    xlabel('h values');
    ylabel('Relative Error');
    title(sprintf('Error Analysis for %s', method_name));
    legend show;
    grid on;
    hold on;
end

val = 1.5;
h = 1;
for i = 1:length(methods)
    Derivative(methods_name{i}, methods{i}, f, f_actual, val, h)
end