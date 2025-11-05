%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(Phi, c, idx, rho, ~, z, u)
   % Linear quadratic problem
   nx = size(z,1);
   c(1:nx) = c(1:nx) - rho * (z - u);
   x = -Phi * c;
   x = x(1:nx);
end