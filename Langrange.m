%% Assignemnt
% Lagrange's Interpolation

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Depatment   : Applied Mechanics


%% Lagrange's Interpolation Function
clc
clear

function [sum] = LI(xs, ys, number)
    % Straight forward brute force way to find l
    l = [];
    
    for j = 1:length(xs)
        a = 1;
        c = 1;
        for i = 1:length(xs)
            if i ~= j       % Or else it will always give 0
                a = a * (number - xs(i));  % Calculating Numerator and Denominator differently
                c = c * (xs(j) - xs(i));
            end
        end
        l(j) = a / c;      
    end
    
    % This calculates l0 x y0 + l1 x y1 + ....
    sum = 0;
    for i = 1: length(l)
        sum = sum + (ys(i) * l(i));
    end
end


function [L] = Lagrange_Basis(xs, j, number)
    L = ones(size(number));
    n = length(xs)
    for i = 1:n
        if i ~= j
            L = L .* (number - xs(i)) / (xs(j) - xs(i))
        end
    end
end

x = [4.0, 5.0, 6.0, 7.0, 8.0]
y = [1.58740105, 1.709976, 1.81712059, 1.912931, 2.0]
sample_points = 50;

% Predicting
test_xs = linspace(min(x), max(x), sample_points);

test_ys = zeros(1, sample_points);
for i = 1:sample_points
    test_ys(i) = LI(x, y, test_xs(i));
end

for j = 1: length(x)
    lb = Lagrange_Basis(x, j, test_xs);
    plot(test_xs, lb, 'LineWidth', 1.5, 'DisplayName', sprintf('P_{%d}(x)', j));
    hold on
end

% Plotting predicted Data
plot(test_xs, test_ys, '--', 'LineWidth', 1.5, 'DisplayName', 'Predicted Data points');
xlabel('X values');
ylabel('Y values');
title("Lagrange's Interpolation", 'FontSize', 25);



hold on

% Plotting actual Data
plot(x, y, 'LineWidth', 1, 'DisplayName', 'Actual Data points', Marker='x', MarkerSize=12);

legend show
grid on;