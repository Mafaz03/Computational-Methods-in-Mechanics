%% Assignemnt
% Newton's Interpolation

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Depatment   : Applied Mechanics



%% Function that finds slope
clc
clear

function [slope] = divided_difference(y2, y1, x2, x1)
    slope = (y2 - y1) / (x2 - x1);

end

%% Newton's Interpolation
function [sum] = NI(x, y, number)

    % The idea:
    % Intead of using matrix to store all the data, we use a single vector
    % and overwrite it, since the non diagonol hold no value to us for this
    % problem, I overwrite the y array itself



    as = [];         % a's are array of the coefficients
    as = [as, y(1)]; % First coefficient is y's first value itself
    
    temp_y = y;      % temp_y is a copy of y, but temp_y keeps shrinking its size, see line 32
    
    for order = 1: length(x)-1       % Number of Columns
        for i = 1: length(temp_y)-1  % Number of Rows
            temp_y(i) = divided_difference(temp_y(i+1), temp_y(i), x(i + order), x(i));   
                                    % Finding Slope, tricky part is the x's where we need to skip ith order of x
        end
        
        temp_y = temp_y(1: end-1); % Shrinking
        as = [as, temp_y(1)];      % Appending to the as array
    end
    
    % This is to compute a0 + a1(x-x0) + a2(x-x0)(x-x1) .....
    sum = 0;
    for i = 1: length(as)
        mul = 1;
        for j = 1:i-1
            mul = mul * (number - x(j));
        end
        sum = sum + (as(i)*mul);
    end
end


%% Using the function
x = [1,  2,    3,    3.2,    3.9];
y = [1,  5,    2,    7,      4];

sample_points = 50;

% Predicting
test_xs = linspace(1, x(end), sample_points);
test_ys = [];
for i = 1: sample_points
    test_ys = [test_ys, NI(x, y, test_xs(i)) ];
end

% Plotting predicted Data
plot(test_xs, test_ys, '--b', 'LineWidth', 1.5, 'DisplayName', 'Predicted Data points');
xlabel('X values');
ylabel('Y values');
title("Newton's Interpolation", 'FontSize', 25);

hold on

% Plotting actual Data
plot(x, y, 'g', 'LineWidth', 1.5, 'DisplayName', 'Actual Data points', Marker='x', MarkerSize=20);

legend show
grid on;