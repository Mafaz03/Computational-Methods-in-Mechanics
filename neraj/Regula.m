clc;
clear;
close all;

f = @(x) sin(10*x) + cos(3*x);
a = 3;
b = 6;

tolerance = 1e-4;
loop = 0;

root_true = 5.67903;

% Ensure the initial bracket is valid
if f(a) * f(b) > 0
    error('f(a) and f(b) must have opposite signs.');
end

% Initial c
c = b - (f(b) * (b - a)) / (f(b) - f(a));

iter_arr  = []; err_arr  = [];

while abs(f(c)) > tolerance
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    
    c = b - (f(b) * (b - a)) / (f(b) - f(a));
    loop = loop + 1;
    iter_arr(loop) = loop;
    err_arr(loop)  = abs(c - root_true);
    fprintf('Loop: %d | c: %.8f | f(c): %.8e | err: %.8e\n', loop, c, f(c), abs(c - root_true));
end

fprintf('\nRoot ≈ %.8f found in %d iterations\n', c, loop);
semilogy(iter_arr, err_arr, 'b-o', 'LineWidth',1.2); hold on;
xlabel('Iteration');
ylabel('True Error |c - root|');
title('True Error vs Iteration');


% --- Estimate order R from the last few iterations ---
% R = [ ln(e_n+1 / e_n) ] / [ ln(e_n / e_n-1) ]

Rvals = zeros(1,length(err_arr)-2);

for k = 3:length(err_arr)
    Rvals(k-2) = log(err_arr(k)/err_arr(k-1)) / log(err_arr(k-1)/err_arr(k-2));
end

R_est = mean(Rvals);  % mean of last few
fprintf('\n\nEstimated order R = %.4f\n', R_est);
