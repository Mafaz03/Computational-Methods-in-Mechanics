%% Assignemnt
% Lagrange Interpolation

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc
clear
close all

%% Part 1 (Preprocessing)
x = [4.0, 5.0, 6.0, 7.0, 8.0];
y = [1.58740105, 1.709976, 1.81712059, 1.912931, 2.0];

n = length(x);
xx = linspace(min(x), max(x), 500);

figure;
hold on;
plot(x, y, 'LineWidth', 1, 'DisplayName', 'Actual Data points', Marker='x', MarkerSize=12);

%% Part 2 (Processing / Using the Interpolation Algorithm)
for order = 1:n-1
    % Take only the first (order+1) data points.
    % For example: 
    % order = 1 -> first 2 points
    % order = 2 -> first 3 points, etc.
    xi = x(1:order+1);
    yi = y(1:order+1);

    yy = zeros(size(xx));
    
    % Loop through each evaluation point xx(k)
    for k = 1:length(xx)
        val = 0;                    % this will accumulate the polynomial value at xx(k)
        
                                    % Build the Lagrange interpolation polynomial term by term
        for i = 1:length(xi)        % Start with basis polynomial L_i(xx(k)) = 1
            L = 1;
            for j = 1:length(xi)
                if j ~= i           % Multiply the factors for L_i(xx(k))
                    L = L * (xx(k) - xi(j)) / (xi(i) - xi(j));
                end
            end
                                    % Add contribution from data point (xi(i), yi(i))
            val = val + yi(i) * L;
        end
        
        yy(k) = val;
    end
    
    plot(xx, yy, 'DisplayName', ['Order ', num2str(order)]);
    hold on;
end

%% Part 3 (Post Processing / Plotting)
% Plotting actual Data
% plot(x, y, 'LineWidth', 1, 'DisplayName', 'Actual Data points', Marker='x', MarkerSize=12);

legend show;
title('Lagrange Interpolation for Orders 0 to n-1');
xlabel('x');
ylabel('y');
grid on;
hold off;
