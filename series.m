%% Assignemnt
% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Depatment   : Applied Mechanics

clc
clear

plot_arr = [];

number = 5;
tolerence = 5e-10;

sum = 1;
loop_completed = 0;
maximum_loops = 400; % So i can break out if the code goes to an infinite Loop, mostly for debugging

actual = exp(number);
i = 1;

while 1
    sum = sum + power(number, i) / factorial(i);

    relative_error = abs((actual - sum)/actual);
    plot_arr = [plot_arr, relative_error];

    if (relative_error) < tolerence;
        % fprintf("\nReached\n")
        break
    end

    loop_completed = loop_completed + 1;

    if loop_completed >= maximum_loops
        break
    end

    i = i + 1;
end

fprintf("relative_error: %.10f\n", relative_error)
fprintf("Loops taken   : %d\n", loop_completed)
fprintf("Predicted     : %.10f\n", sum)
fprintf("Actual        : %.10f\n", exp(number))

% length(plot_arr)
plot(1:length(plot_arr), plot_arr);
title('Error vs Itterations')
xlabel("Itteration")
ylabel("Realtive erorr")
