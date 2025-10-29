% Student Details
% Roll Number : AM25M017
% Name        : Y.Sujith
% Assignment  : 8
% 2D HEAT CONDUCTION WITH HEAT GENERATION Q1

clear all; clc; close all;

Width = 1;
Height = 1.5;
tolerance = 1e-6;
maxIterations = 10000;
thermalConductivity = 0.4;
heatSource = 100 / thermalConductivity;
Nx_values = [9 11 20];

% LIEBMANN 5 POINTS (SOR)
for gridIdx = 1:length(Nx_values)
    Nx = Nx_values(gridIdx);
    deltaH = Width / (Nx - 1);
    Ny = round(Height / deltaH + 1);
    
    rootValue = cos(pi / (Nx - 1)) + cos(pi / (Ny - 1));
    omegaSOR = 4 / (2 + sqrt(4 - rootValue^2));
    
    TempField = zeros(Ny, Nx);
    TempField(1, :) = 0;
    TempField(Ny, :) = 0;
    
    tic;
    for iteration = 1:maxIterations
        maxError = 0;
        
        for row = 2 : Ny - 1
            for col = 1 : Nx
                oldTemp = TempField(row, col);
                
                if col == 1
                    newTemp = 0.25 * (2 * TempField(row, 2) + TempField(row+1, col) + TempField(row-1, col) + heatSource * deltaH^2);
                elseif col == Nx
                    newTemp = 0.25 * (2 * TempField(row, Nx-1) + TempField(row+1, col) + TempField(row-1, col) + heatSource * deltaH^2);
                else
                    newTemp = 0.25 * (TempField(row, col+1) + TempField(row, col-1) + TempField(row+1, col) + TempField(row-1, col) + heatSource * deltaH^2);
                end
                
                TempField(row, col) = oldTemp + omegaSOR * (newTemp - oldTemp);
                
                errorVal = abs(oldTemp - TempField(row, col));
                if errorVal > maxError
                    maxError = errorVal;
                end
            end
        end
        
        if maxError < tolerance
            break;
        end
    end
    timer = toc;
    
    fprintf('LB 5-Point For h = %.5f, %d iterations \n %.4f seconds \n Error = %5e \n', deltaH, iteration, timer, maxError);
    
    xAxis = linspace(0, Width, Nx);
    yAxis = linspace(0, Height, Ny);
    figure(gridIdx);
    subplot(1, 2, 1);
    surf(xAxis, yAxis, TempField);
    title('LB 5 Points Final Temperature (Surface)');
    xlabel('x'); ylabel('y');
    subplot(1, 2, 2);
    contourf(xAxis, yAxis, TempField, 20);
    title('LB 5 Points Final Temperature (Contour)');
    xlabel('x'); ylabel('y');
end

%---------------------------------- LIEBMANN (SOR) 9 POINTS -------------------------------------------
disp('-----------------------------------------LIEBMANN 9 POINTS (SOR)-------------------------------------------------------------')

for gridIdx = 1:length(Nx_values)
    Nx = Nx_values(gridIdx);
    deltaH = Width / (Nx - 1);
    Ny = round(Height / deltaH + 1);
    
    rootValue = cos(pi / (Nx - 1)) + cos(pi / (Ny - 1));
    omegaSOR = 4 / (2 + sqrt(4 - rootValue^2));
    
    TempField = zeros(Ny, Nx);
    TempField(1, :) = 0;
    TempField(Ny, :) = 0;
    
    tic;
    for iteration = 1:maxIterations
        maxError = 0;
        
        for row = 2 : Ny - 1
            for col = 1 : Nx
                oldTemp = TempField(row, col);
                
                if col == 1
                    newTemp = 0.05 * ( 4*(2*TempField(row, 2) + TempField(row+1, col) + TempField(row-1, col)) ...
                        + (2*TempField(row+1, 2) + 2*TempField(row-1, 2)) + 6*heatSource*deltaH^2 );
                elseif col == Nx
                    newTemp = 0.05 * ( 4*(2*TempField(row, Nx-1) + TempField(row+1, col) + TempField(row-1, col)) ...
                        + (2*TempField(row+1, Nx-1) + 2*TempField(row-1, Nx-1)) + 6*heatSource*deltaH^2);
                else
                    newTemp = 0.05 * ( 4*(TempField(row, col+1) + TempField(row, col-1) + TempField(row+1, col) + TempField(row-1, col)) ...
                        + (TempField(row+1, col+1) + TempField(row-1, col+1) + TempField(row+1, col-1) + TempField(row-1, col-1)) ...
                        + 6*heatSource*deltaH^2 );
                end
                
                TempField(row, col) = oldTemp + omegaSOR * (newTemp - oldTemp);
                
                errorVal = abs(oldTemp - TempField(row, col));
                if errorVal > maxError
                    maxError = errorVal;
                end
            end
        end
        if maxError < tolerance
            break;
        end
    end
    timer = toc;
    
    fprintf('LB 9-Point For h = %.5f, %d iterations \n %.4f seconds \n Error = %5e \n', deltaH, iteration, timer, maxError);
    
    xAxis = linspace(0, Width, Nx);
    yAxis = linspace(0, Height, Ny);
    figure(gridIdx + 3);
    subplot(1, 2, 1);
    surf(xAxis, yAxis, TempField);
    title('LB 9 Points Final Temperature (Surface)');
    xlabel('x'); ylabel('y');
    subplot(1, 2, 2);
    contourf(xAxis, yAxis, TempField, 20);
    title('LB 9 Points Final Temperature (Contour)');
    xlabel('x'); ylabel('y');
end

%---------------------------------- GAUSS-SEIDEL 5 POINTS -------------------------------------------
disp('-----------------------------------GAUSS SIEDEL 5 POINTS-------------------------------------------------------------------')

for gridIdx = 1:length(Nx_values)
    Nx = Nx_values(gridIdx);
    deltaH = Width / (Nx - 1);
    Ny = round(Height / deltaH + 1);
    
    TempGS5 = zeros(Ny, Nx);
    TempGS5(1, :) = 0;
    TempGS5(Ny, :) = 0;
    
    tic;
    for iteration = 1:maxIterations
        maxError = 0;
        for row = 2 : Ny - 1
            for col = 1 : Nx
                oldTemp = TempGS5(row, col);
                
                if col == 1
                    TempGS5(row, col) = 0.25 * (2*TempGS5(row, 2) + TempGS5(row+1, col) + TempGS5(row-1, col) + heatSource*deltaH^2);
                elseif col == Nx
                    TempGS5(row, col) = 0.25 * (2*TempGS5(row, Nx-1) + TempGS5(row+1, col) + TempGS5(row-1, col) + heatSource*deltaH^2);
                else
                    TempGS5(row, col) = 0.25 * (TempGS5(row, col+1) + TempGS5(row, col-1) + TempGS5(row+1, col) + TempGS5(row-1, col) + heatSource*deltaH^2);
                end
                
                errorVal = abs(oldTemp - TempGS5(row, col));
                if errorVal > maxError
                    maxError = errorVal;
                end
            end
        end
        if maxError < tolerance
            break;
        end
    end
    
    disp(TempGS5)
    timer = toc;
    
    fprintf('GS 5-Point For h = %.5f, %d iterations \n %.4f seconds \n Error = %5e \n', deltaH, iteration, timer, maxError);
    
    xAxis = linspace(0, Width, Nx);
    yAxis = linspace(0, Height, Ny);
    figure(gridIdx + 6);
    subplot(1, 2, 1);
    surf(xAxis, yAxis, TempGS5);
    title('GS 5 Points Final Temperature (Surface)');
    xlabel('x'); ylabel('y');
    subplot(1, 2, 2);
    contourf(xAxis, yAxis, TempGS5, 20);
    title('GS 5 Points Final Temperature (Contour)');
    xlabel('x'); ylabel('y');
end

%---------------------------------- GAUSS-SEIDEL 9 POINTS -------------------------------------------
disp('-----------------------------------------GAUSS SIEDEL 9 POINTS-------------------------------------------------------------')

for gridIdx = 1:length(Nx_values)
    Nx = Nx_values(gridIdx);
    deltaH = Width / (Nx - 1);
    Ny = round(Height / deltaH + 1);
    
    TempGS9 = zeros(Ny, Nx);
    TempGS9(1, :) = 0;
    TempGS9(Ny, :) = 0;
    
    tic;
    for iteration = 1:maxIterations
        maxError = 0;
        for row = 2 : Ny - 1
            for col = 1 : Nx
                oldTemp = TempGS9(row, col);
                
                if col == 1
                    TempGS9(row, col) = 0.05 * (4*(2*TempGS9(row, 2) + TempGS9(row+1, col) + TempGS9(row-1, col)) ...
                        + (2*TempGS9(row+1, 2) + 2*TempGS9(row-1, 2)) + 6*heatSource*deltaH^2);
                elseif col == Nx
                    TempGS9(row, col) = 0.05 * (4*(2*TempGS9(row, Nx-1) + TempGS9(row+1, col) + TempGS9(row-1, col)) ...
                        + (2*TempGS9(row+1, Nx-1) + 2*TempGS9(row-1, Nx-1)) + 6*heatSource*deltaH^2);
                else
                    TempGS9(row, col) = 0.05 * (4*(TempGS9(row, col+1) + TempGS9(row, col-1) + TempGS9(row+1, col) + TempGS9(row-1, col)) ...
                        + (TempGS9(row+1, col+1) + TempGS9(row-1, col+1) + TempGS9(row+1, col-1) + TempGS9(row-1, col-1)) ...
                        + 6*heatSource*deltaH^2);
                end
                
                errorVal = abs(oldTemp - TempGS9(row, col));
                if errorVal > maxError
                    maxError = errorVal;
                end
            end
        end
        if maxError < tolerance
            break;
        end
    end
    
    disp(TempGS9)
    timer = toc;
    
    fprintf('GS 9-Point For h = %.5f, %d iterations \n %.4f seconds \n Error = %5e \n', deltaH, iteration, timer, maxError);
    
    xAxis = linspace(0, Width, Nx);
    yAxis = linspace(0, Height, Ny);
    figure(gridIdx + 9);
    subplot(1, 2, 1);
    surf(xAxis, yAxis, TempGS9);
    title('GS 9 Points Final Temperature (Surface)');
    xlabel('x'); ylabel('y');
    subplot(1, 2, 2);
    contourf(xAxis, yAxis, TempGS9, 20);
    title('GS 9 Points Final Temperature (Contour)');
    xlabel('x'); ylabel('y');
end
