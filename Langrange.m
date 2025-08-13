clc
clear

xs = [0, 0.6, 1.2];      
ys = [1, 3, 2];          

x = 0.7;
l = [];

for j = 1:length(xs)
    a = 1;
    c = 1;
    for i = 1:length(xs)
        if i ~= j
            a = a * (x - xs(i));
            c = c * (xs(j) - xs(i));
        end
    end
    l(j) = a / c; 
end

sum = 0;
for i = 1: length(l)
    sum = sum + (ys(i) * l(i));
end


disp(sum)