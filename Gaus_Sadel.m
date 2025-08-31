%% Assignemnt
% Gaus Sadel Itterator

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;

loops_taken = 0;

A = [4   1  -1;     % A Matrix
     1  -5  -1; 
     2  -1  -6];

B = [13; -8; -2];   % B Matrix

X = [0; 0; 0];      % Initial Guess
X_new = [0; 0; 0];

tolerence = 1e-12;
relative_error = 0;

%% Part 1 (Preprocessing)
% Checking if Diagonally Dominant
function [] = Diag_dom(A)
    diag_dom = 0;
    for j = 1:length(A)
        sum = 0;
        for i = 1:length(A)
            if i ~= j
                sum = sum + abs(A(j, i)); % abs value of sum of non diagnoal elements
            end
        end
        if sum > abs(A(j, j))             % Comparing to diagnol elements
            diag_dom = diag_dom + 1
        end
    end
    
    if diag_dom > 0
        fprintf('Matrix is not Diagonally dominant :(\n\n')
    else
        fprintf('Matrix is Diagonally dominant :)\n\n')
    end
end

Diag_dom(A)

%% Part 2 (Processing)
while (relative_error > tolerence) || loops_taken == 0
    for j = 1:length(A)
        sum = 0;
        for i = 1:length(A)
            if i ~= j
                sum = sum + (A(j, i) * X_new(i));            % Sum of non diagnol elements
            end
        end
        X_new(j) = (B(j) - sum) / A(j, j);                   % b - sum
                                                             % Gauss sadel uses new values right away in the same loop
    end

    relative_error = max(abs(X_new - X) ./ (X_new + 1e-9));  % Calculating relative error

    X = X_new;                                               % Copying new array to actual array of initial guesses
    loops_taken = loops_taken + 1;
    
    %% Part 3 (post processing)
    % Printing it out
    fprintf("Loop: %d   |", loops_taken)
    for i = 1:length(X)
        fprintf("   X_%d: %f   |   ", i, X(i))
    end
    fprintf("\n")
end

fprintf("\nLoops Taken: %d\n", loops_taken)
fprintf("\nRelative Error: %d\n", relative_error)

