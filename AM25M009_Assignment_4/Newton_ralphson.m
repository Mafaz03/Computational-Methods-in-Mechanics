%% Assignemnt - 4
% Newton Ralphson

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;
close all;

%% Part 1 (Preprocessing)
x0 = 10;
tol = 1e-6;
max_loop = 100;

root_true = 1;

syms x
f_sym = (x-1) * (exp(x-1) - 1);
% f_sym = (x-1)*(x-2);
df_sym = diff(f_sym);

% Convert back to function 
f = matlabFunction(f_sym);
df = matlabFunction(df_sym);

rel_error = inf;
loop = 0;

iter_arr = [];  err_arr = [];

%% Part 2 (Processing)
while rel_error > tol
    x1 = x0 - f(x0)/df(x0);
    
    rel_error = abs((x1 - x0) / x1);
    
    x0 = x1;
    loop = loop + 1;

    iter_arr(loop) = loop;
    err_arr(loop)  = abs(x1 - root_true);

   
    fprintf("Loop: %d | c: %f\n", loop, x1)
end


%% Part 3 (Post Processing / Plotting)
semilogy(iter_arr, err_arr, 'b-o', 'LineWidth',1.2); hold on;
xlabel('Iteration');
ylabel('True Error |c - root|');
title('True Error vs Iteration');


% R = [ ln(e_n+1 / e_n) ] / [ ln(e_n / e_n-1) ]

Rvals = zeros(1,length(err_arr)-2);

for k = 3:length(err_arr)
    % Rvals(k-2) = log(err_arr(k)/err_arr(k-1)) / log(err_arr(k-1)/err_arr(k-2));
    Rvals(k-2) = abs(log(err_arr(k)/err_arr(k-1)) / log(err_arr(k-1)/err_arr(k)));
end

R_est = mean(Rvals);  % mean of last few
fprintf('\n\nEstimated order R = %.4f\n', R_est);

