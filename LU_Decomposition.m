%% Assignment
% LU Decomposition (Doolittle method)

% Name        : Mohamed Mafaz
% Roll Number : AM25M009
% Department  : Applied Mechanics

clc;
clear;

%% Part 1 (Preprocessing)
A = [2  -1  3   2;
     2   2  0   4;
     1   1 -2   2;
     1   3  4  -1];



%% Part 2 (Implementing)
function [L, U] = LU(A)
    n = size(A,1);
    L = eye(n);
    U = A;
    for i = 1:n-1
        for j = i+1:n
            m = U(j,i)/U(i,i);           % multiplier
            L(j,i) = m;                  % store multiplier in L
            U(j,:) = U(j,:) - m*U(i,:);  % eliminate in U
        end
    end
end

%% Part 3 (Post Processing)
[L, U] = LU(A)
fprintf('L matrix:\n')
disp(L)

fprintf('U matrix:\n')
disp(U)

fprintf('A matrix:\n')
disp(round(L * U, 1))