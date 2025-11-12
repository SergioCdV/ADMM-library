%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(Fv, V, c, rho, ~, z, u)
   % Linear quadratic problem
   res = rho * (z - u);
   x = Fv .* ( c - V * res );
end