%% Student Details
% Roll number : AM25M009
% Name        : Mohamed Mafaz
% Assignment  : Maclaurin Series and error approximation
% Department  : Applied Mechanics


clc;
clear;

number = 0.2*pi;    % x: where we want to find e(x)


plot_arr = [];      % To plot relative error wrt to itterations

tolerence = 5e-9;

sum = 1;
loop_completed = 0;
maximum_loops = 100; % So i can break out if the code goes to an infinite Loop, mostly for debugging

actual = exp(number);


%% Part 1 (Preprocessing / Writing Serie's loop)
% Maclaurin Series
i = 1;     % Starting from 2nd factor, since 1 is always present

while 1
    sum = sum + power(number, i) / factorial(i);  % maclaurin series 

    relative_error = abs((actual - sum)/actual);  % Relative error: |x_true - x| / x_true
    plot_arr = [plot_arr, relative_error];        % Storing error since i want to plot it

    if (relative_error) < tolerence;              % Comparing float is generally not a good idea for extremly small numbers due to machine precision 
                                                  % But since error arent that small it wont cause issues.
        break
    end

    loop_completed = loop_completed + 1;          % Loops are kept track of

    if loop_completed >= maximum_loops            % Refer line 16
        break
    end

    i = i + 1;
end

%% Part 2 (Processing / Using the Loop)
% Printing Stats
fprintf("relative_error: %.10f\n", relative_error)
fprintf("Terms used   : %d\n", loop_completed+1)
fprintf("Predicted     : %.10f\n", sum)
fprintf("Actual        : %.10f\n", exp(number))


%% Part 3 (post processing or plots or results)
% Plotting Relative error
plot(1:length(plot_arr), plot_arr, '-o');
title('Error vs Itterations', 'FontSize', 25)
xlabel("Number of Terms")
ylabel("Realtive erorr")
