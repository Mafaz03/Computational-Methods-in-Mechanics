%% Student Details
clc;
clear;

A = [10 -1 2; 
    -1 11 -1; 
    2 -1 10];

B = [6 25 -11];

X = [0 0 0];
X_new = [0 0 0];

tolerence = 1e-60;
relative_error = tolerence + 1;

while relative_error > tolerence
    for j = 1:length(B)
        for i = 1:length(B)
            sum = 0;
            if i ~= j
                sum = sum + (A(j, i) * X_new(i));
            end
        
        end
        X_new(j) = (B(j) - sum) / A(j, j);
    end
    X = X_new;
    relative_error = max(abs(X_new - X) ./ abs(X_new));
end

X
