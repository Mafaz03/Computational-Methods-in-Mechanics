clc;
clear;

f =@(x) exp(-x^2);

a = 0;
b = 1;

%% Single Step

% Trapezoid
h = (b-a);
trap_integral = h/2 * (f(a) + f(b));
fprintf("Trapezoid Single step: %f\n", trap_integral)

% Simpsion's 1/3rd
h = (a + b)/2;
simp_1_3_integral = (h/3) * (f(a) + (4 * f(a + h)) + f(b));
fprintf("Simpsion's 1/3rd Single step: %f\n", simp_1_3_integral)

% Simpsion's 3/8rd
h = (a + b)/3;
simp_3_8_integral = ((3*h)/8) * (f(a) + (3 * f(a + h)) + (3 * f(a + 2*h)) + f(b));
fprintf("Simpsion's 3/8th Single step: %f\n", simp_3_8_integral)

%% Multiple Steps

% Trapezoid

M = 100;
h = (b-a)/M;
k_val = linspace(a, b, M);

sum = 0;
for k = 1:M-1
    sum = sum + (   f(k_val(k))  +  f(k_val(k+1))   );
end

fprintf("\n\nTrapezoid Multi step: %f\n", (h/2) * sum)


% Simpsion's 1/3rd
M = 100;

n = 2 * M; % n = 2M
k_val = linspace(a, b, n);
h = (b-a)/M;
sum = 0;
for k = 1:n
    if mod(k, 2) == 0
        sum = sum + (2 * f(k_val(k)));
    else
        sum = sum + (4 * f(k_val(k)));
    end
end
simp_1_3_integral = (h/3) * (f(k_val(1)) + sum + f(k_val(end)));
fprintf("Simpsion's 1/3rd Multi step: %f\n", simp_1_3_integral)
