%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/12/23
% File: Test_11.m 
% Issue: 0 
% Validated: 

%% Test exercise XI %% 
% Solve for the minimum of 0.5 * x.' * P * x, Ax \in C % 

%% Input data
P = 1e-3 * eye(6);  
P = 0.5 * (P + P.') + 1e-2 * eye(6); 
q = rand(6,1);

A = rand(3,6);

%% Solve the problem 
Solver = Solvers.AQP_solver(@(x)myConeCheck(x), @(x)myConeProjec(x), P, q, A);

[x, z, Output] = Solver.solver();

%% Auxiliary functions 
function [x_proj] = myConeProjec(x)
    x_proj = x;
    for i = 1:size(x)
        if abs(x(i)) > 1
            x_proj(i) = 1 * sign(x(i));
        end
    end
end

function [index] = myConeCheck(x)
    index = x;
    for i = 1:size(x)
        if abs(x(i)) > 1
            index(i) = true;
        else
            index(i) = false;
        end
    end
end




