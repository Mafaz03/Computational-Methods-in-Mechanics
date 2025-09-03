clc
clear

loop = 0;

f = @(x) (x-1) * (x-2) * (x-3);
a = -100;
b = 100;

tolerance = 1e-5;

c = b - ((f(b) * (b - a)) / (f(b) - f(a)))

while abs(f(c)) > tolerance
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
    loop = loop + 1;
    fprintf("Loop: %d | c: %f\n", loop, c)
end
