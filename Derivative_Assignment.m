
%% AM25M009 - Assignment 5 - Theory Class

clc; clear;
close all;

f  = @(x) exp(-2*x) - x;
x0 = 2;

% exact derivative
fprime_exact = -2*exp(-2*x0) - 1;

% step sizes
h_values = 0.5:-0.001:0.05;

fdash_values = zeros(size(h_values));
err_values   = zeros(size(h_values));
fprintf("NOT PRINTING ALL!!\n\n")
fprintf('%8s %15s %15s\n','h','Derivative','Abs_Error');
fprintf('-----------------------------------------------\n');

for k = 1:length(h_values)
    h = h_values(k);
    fdash = ( f(x0 + h) - f(x0 - h) ) / (2*h);     % central difference
    err   = abs(fdash - fprime_exact);

    fdash_values(k) = fdash;
    err_values(k)   = err;

    if mod(k,15) == 0
        fprintf('%8.3f %15.8f %15.8e\n', h, fdash, err);
    end
end

% linear plot
figure;
plot(h_values, err_values, 'o-','LineWidth',1.5);
xlabel('Step size h');
ylabel('Absolute error');
title('Central difference derivative: Error vs h');

% log–log
figure;
loglog(h_values, err_values, 'o-','LineWidth',1.5);
xlabel('Step size h (log scale)');
ylabel('Absolute error (log scale)');
title('Central difference derivative: Log–log error plot');