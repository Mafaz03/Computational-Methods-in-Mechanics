clc; clear;

f  = @(x) exp(-2*x) - x;
x0 = 2;

% exact derivative
fprime_exact = -2*exp(-2*x0) - 1;

% step sizes
h_values = 0.5:-0.001:0.05;

fprintf('%8s %15s %15s\n','h','Derivative','Abs_Error');
fprintf('-----------------------------------------------\n');

for h = h_values
    fdash = ( f(x0 + h) - f(x0 - h) ) / (2*h);
    err   = abs(fdash - fprime_exact);
    fprintf('%8.3f %15.8f %15.8e\n', h, fdash, err);
end