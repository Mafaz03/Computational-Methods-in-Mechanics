f = @(x) x^3 - x - 2;

a = 1;
b = 2;

if f(a) * f(b) > 0 % Condition checking
    error("Make sure a and b give out f(a) and f(b) in the opposite quadrants")