clc
clear

as = [];

function [slope] = divided_difference(y2, y1, x2, x1)
    slope = (y2 - y1) / (x2 - x1);

end


function [sum] = NI(x, y, number)

    as = []

    as = [as, y(1)]
    
    temp_y = y;
    
    for order = 1: length(x)-1
        for i = 1: length(temp_y)-1
            temp_y(i) = divided_difference(temp_y(i+1), temp_y(i), x(i + order), x(i));
        
        end
        
        temp_y = temp_y(1: end-1);
        as = [as, temp_y(1)];
    end
    
    
    sum = 0
    
    for i = 1: length(as)
        mul = 1
        for j = 1:i-1
            mul = mul * (number - x(j));
        end
        sum = sum + (as(i)*mul);
    end
end


x = [1,  2,    3,    3.2,    3.9];
y = [1,  5,    2,    7,      4];

sample_points = 50
test_xs = linspace(1, length(x), sample_points)
test_ys = []
for i = 1: sample_points
    test_ys = [test_ys, NI(x, y, test_xs(i)) ]
end


plot(1:length(test_ys), test_ys, '--b', 'LineWidth', 1.5);
xlabel('X values');
ylabel('Y values');
title('Newtons Interpolation');
legend('Predicted Data points');

hold on

plot(x, y, '--g', 'LineWidth', 1.5);
legend('Actual Data points');


grid on;