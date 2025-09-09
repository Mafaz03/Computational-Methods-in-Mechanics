clc; clear; close all;

% Given functions
f = @(x) exp(-2*x) - x;
fp_actual = @(x) -2*exp(-2*x) - 1;
x = 2;

hs = 10:-0.001:0.05; % Step sizes from 0.15 to 0.05
num_steps = length(hs);
error_central = zeros(1,num_steps);
centralVal = zeros(1,num_steps);

for i = 1:num_steps
    h = hs(i);
    % Central difference calculation
    centralVal(i) = (f(x+h) - f(x-h)) / (2*h);
    error_central(i) = abs(centralVal(i) - fp_actual(x));
end

% Display errors for h=0.5
disp(['Central diff at h=0.5: ', num2str(centralVal(1)), ...
      ', abs error: ', num2str(error_central(1))]);

% Plot error vs h
figure;
plot(hs, error_central, 'LineWidth',2);
xlabel('Step size h');
ylabel('Absolute error');
title('Central difference absolute error vs h');
grid on;

% Find the optimal step size (where error is minimized)
[~, opt_idx] = min(error_central);
optimal_h = hs(opt_idx);
fprintf('Optimal step size is approximately h = %.4f with minimum abs error = %e\n', ...
        optimal_h, error_central(opt_idx));
