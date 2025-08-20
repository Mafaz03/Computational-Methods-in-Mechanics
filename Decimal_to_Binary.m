% Student Details
% Roll number : AM25M009
% Name        : Mohamed Mafaz
% Assignment  : Decimal to Binary up to 8-bit precision
% Department  : Applied Mechanics


clc;
clear;

%% Part 1 (Preprocessing / Writing Functions)

function [] = D2B(number, precision)
    
    before_decimal = floor(number);
    after_decimal  = number - before_decimal;   % Splitting the Decimal into 2 parts -> before_decimal, after_decimal
    before_decimal_list = [];

    % Doing Calculations to digits before decimal seperately
    while before_decimal ~= 0
        remainder        = mod(before_decimal, 2);       % 0 or 1
        before_decimal   = floor(before_decimal / 2);    % before_decimal is not the Quotient
        before_decimal_list   = [before_decimal_list, remainder];  % Adding to the array
    
    end
    
    decimal_after_list = [];
    
    % Doing Calculations to digits after decimal seperately
    digit_count = 0;
    if after_decimal ~= 0
        while after_decimal
            a = after_decimal * 2;                   % Multiply with 2 and see if the digit before the decimal is 0 or 1
            binary = floor(a);
            decimal_after_list = [decimal_after_list, binary]; % Adding 0 or 1 to the list
            after_decimal = a - binary;              % Subtracting a with 0 or 1

            digit_count = digit_count + 1;
            if digit_count >= precision             % Incase of never ending loop, we can exit       
                break
            end
        end 
          
    end


    fprintf("Decimal: %f | Binary: ", number)

    % Printing the before_decimal_list in reverse order, since that is the algorithm
    for i = 0:length(before_decimal_list)-1
        fprintf("%d", before_decimal_list(end-i));
    end
    % if binary starts after decimal point, then instead of .101 we will use 0.101
    if length(before_decimal_list) == 0
        fprintf("0")
    end
    % Only print '.' if there is binary after the decimal point also
    if length(decimal_after_list) ~= 0
        fprintf(".")
    end
    % Print after decimal point binary characters if there are any
    for i = 1:length(decimal_after_list)
        fprintf("%d", decimal_after_list(i));
    end

    if digit_count >= precision
        fprintf(" (Precision reached: %d digits)", precision);
    end
    
    fprintf("\n")
end


%% Part 2 (Processing / Using the function)

numbers = [5.625, 0.8925, 205, 124.456];

for num = numbers
    if isempty(num) || ~isnumeric(num)  % Only Execute if user enters a number
        disp('No number entered. Exiting...');
    
    else
        D2B(num, 8)
    end
end