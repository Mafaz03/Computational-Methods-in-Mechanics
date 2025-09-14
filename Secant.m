clc;
clear;
close all;

x0 = 10;
x1 = 9; 
tol = 1e-6;

root_true = 1;

syms x
f = @(x) (x-1) * (exp(x-1) - 1);

rel_error = inf;
loop = 0;
iter_arr = [];
err_arr = [];

while rel_error > tol
    x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
    
    rel_error = abs((x2 - x1) / x2);
    
    % Update values for next iteration
    x0 = x1;
    x1 = x2;
    loop = loop + 1;
    
    iter_arr(loop) = loop;
    err_arr(loop) = abs(x2 - root_true);
    
    fprintf("Loop: %d | c: %f\n", loop, x2);
end

% Plot true error in semilog scale
semilogy(iter_arr, err_arr, 'b-s', 'LineWidth',1.2);
xlabel('Iteration');
ylabel('True Error |c - root|');
title('True Error vs Iteration');

% Estimate order R from last few iterations
Rvals = zeros(1,length(err_arr)-2);
for k = 3:length(err_arr)
    Rvals(k-2) = log(err_arr(k)/err_arr(k-1)) / log(err_arr(k-1)/err_arr(k-2));
end

R_est = mean(Rvals);  % mean of last few
fprintf('\n\nEstimated order R = %.4f\n', R_est);
