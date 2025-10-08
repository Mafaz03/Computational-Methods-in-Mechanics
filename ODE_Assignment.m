clc;
clear;

% ODE definition and exact solution
f_dash = @(t, y) -1000 * (y - (t + 2)) + 1;
exact_sol = @(t) -exp(-1000*t) + t + 2;

% Step sizes to try
h_values = [5e-4, 2e-3, 2.5e-3];
t_final = 0.01;

for k = 1:length(h_values)
    h = h_values(k);
    t = 0:h:t_final;
    y = zeros(size(t));
    y(1) = 1; % initial condition
    
    err = zeros(size(t));
    y_exact = zeros(size(t));
    for i = 1:length(t)-1
        y(i+1) = y(i) + h * f_dash(t(i), y(i));
        y_exact(i+1) = exact_sol(t(i+1));
        err(i+1) = abs(y(i+1) - y_exact(i+1));
    end

    
    % Plot solution and error for each h
    figure(k)
    subplot(2,1,1)
    plot(t,y,'b.-',t,y_exact,'r--');
    legend('Euler','Exact')
    xlabel('t'), ylabel('y')
    title(['Euler method vs. exact, h = ', num2str(h)])
    grid on

    fprintf("Step size h = %.4g | y_Euler(0.01) = %f | y_exact(0.01) = %f | error at t=0.01 = %e\n", ...
        h, y(end), y_exact(end), err(end));
    
end
