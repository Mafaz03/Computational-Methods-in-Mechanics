%% Assignemnt - 4
% Bisection

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc
clear
close all

%% Part 1 (Preprocessing)

loop = 0;
f = @(x) sin(10*x) + cos(3*x);
a = 3;
b = 6;

tolerance = 1e-4;
root_true = 3.74575;

iter_arr = [];  err_arr = [];

%% Part 2 (Processing)
c = (b + a)/2;
while abs(f(c)) > tolerance
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    
    c = (b + a) / 2;
    loop = loop + 1;
    iter_arr(loop) = loop;
    err_arr(loop)  = abs(c - root_true);

    fprintf("Loop: %d | c: %f\n", loop, c)
end

semilogy(iter_arr, err_arr, 'b-o', 'LineWidth',1.2); hold on;
xlabel('Iteration');
ylabel('True Error |c - root|');
title('True Error vs Iteration');

%% Part 3 (Post Processing / Plotting)

% R = [ ln(e_n+1 / e_n) ] / [ ln(e_n / e_n-1) ]

Rvals = zeros(1,length(err_arr)-2);

for k = 3:length(err_arr)
    Rvals(k-2) = log(err_arr(k)/err_arr(k-1)) / log(err_arr(k-1)/err_arr(k-2));
end

R_est = mean(Rvals);  % mean of last few
fprintf('\n\nEstimated order R = %.4f\n', R_est);


