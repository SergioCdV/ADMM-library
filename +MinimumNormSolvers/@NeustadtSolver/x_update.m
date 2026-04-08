%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(Phi, A, q, b, rho, ~, z, u)
    % Linear quadratic problem
    v = z - u;
    
    opts.LT = true;
    lambda  = rho * (A * v - b) - A * q;
    lambda  = linsolve(Phi, lambda, opts);
    
    opts.LT = false; 
    opts.UT = true;
    lambda  = linsolve(Phi.', lambda, opts);
    
    x = v - (q + A.' * lambda) / rho;
end
