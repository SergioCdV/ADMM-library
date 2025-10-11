%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Proximal minimization for X
function [x] = x_update(Phi, b, rho, x, z, u)
    % Impulses update
    e = [rho * (z - u); b]; 
    sol = Phi * e;

    x = sol(1:length(x));
end