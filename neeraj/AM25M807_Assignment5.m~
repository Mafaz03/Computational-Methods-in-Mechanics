%Roll number : AM25M807
%Name : NIRAJ KUMAR SINGH
%Assignment: Numerical Derivative Assighnment 5 (Q1) 
% -------------------------------
%Clear workspace and command window


clc;
clear all;
close all;

% Define symbolic variable and test function
syms t
test_func = exp(t) * cos(t);

% Create function handles for numerical evaluation
func_handle = matlabFunction(test_func);
derivative_exact = matlabFunction(diff(test_func, t));

% Define numerical differentiation schemes
diff_forward  = @(func, point, step) (func(point + step) - func(point)) / step;
diff_backward = @(func, point, step) (func(point) - func(point - step)) / step;
diff_centered = @(func, point, step) (func(point + step) - func(point - step)) / (2 * step);

% Store methods in cell array with corresponding names
numerical_methods = {diff_forward, diff_backward, diff_centered};
scheme_labels = ["Forward Difference", "Backward Difference", "Central Difference"];

% Analysis function for each differentiation method
function [] = AnalyzeMethod(label, diff_scheme, target_func, exact_deriv, eval_point, initial_step)
    step_sizes = [];
    error_values = [];
    
    fprintf("\n\n======== %s Analysis ========\n\n", label)
    
    error_threshold = 1e-6;
    iteration = 0;
    current_error = inf;
    current_step = initial_step;
    
    while (current_error > error_threshold) || iteration == 0
        % Calculate numerical derivative
        numerical_result = diff_scheme(target_func, eval_point, current_step);
        exact_result = exact_deriv(eval_point);
        
        % Compute relative error
        current_error = abs((numerical_result - exact_result) / exact_result);
        
        % Store data for plotting
        step_sizes = [step_sizes, current_step];
        error_values = [error_values, current_error];
        
        % Update step size and iteration counter
        current_step = current_step / 1.5;
        iteration = iteration + 1;
        
        fprintf("%d | %.8f | %.8f | error = %.2e\n", iteration, numerical_result, exact_result, current_error);
    end
    
    % Create error plot
    plot(error_values, step_sizes)
    xlabel('Step Size (h)');
    ylabel('Relative Error');
    title(sprintf('Convergence Analysis: %s', label));
    legend show;
    grid on;
    hold on;
end

% Set evaluation parameters
evaluation_point = 1.5;
starting_step = 1;

% Execute analysis for each method
for method_idx = 1:length(numerical_methods)
    AnalyzeMethod(scheme_labels{method_idx}, numerical_methods{method_idx}, ...
                  func_handle, derivative_exact, evaluation_point, starting_step)
end