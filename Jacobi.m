%% Student Details
clc;
clear;

loops_taken = 0;

A = [4 -1 0; 
    -1 4 -1; 
    0 -1 3];

B = [15; 10; 10];

X = [0; 0; 0];
X_new = [0; 0; 0];

tolerence = 1e-9;
relative_error = 0;

% Checking if Diagonally Dominant
function [] = Diag_dom(A)
    diag_dom = 0;
    for j = 1:length(A)
        for i = 1:length(A)
            sum = 0;
            if i ~= j
                sum = sum + (A(j, i));
            end
        end
        if sum > A(j, j)
            diag_dom = diag_dom + 1
        end
    end
    
    if diag_dom > 0
        fprintf('Matrix is not Diagonally dominant.\n\n')
    else
        fprintf('Matrix is Diagonally dominant :)\n\n')
    end
end

Diag_dom(A)

while (relative_error > tolerence) || loops_taken == 0
    for j = 1:length(A)
        sum = 0;
        for i = 1:length(A)
            if i ~= j
                sum = sum + (A(j, i) * X_new(i));
            end
        end
        X_new(j) = (B(j) - sum) / A(j, j);
    end

    relative_error = max(abs(X_new - X) ./ (X_new + 1e-9));

    X = X_new;
    loops_taken = loops_taken + 1;
end


for i = 1:length(X)
    fprintf("X_%d: %f\n", i, X(i))
end
fprintf("\nLoops Taken: %d\n", loops_taken)
fprintf("\nRelative Error: %d\n", relative_error)

