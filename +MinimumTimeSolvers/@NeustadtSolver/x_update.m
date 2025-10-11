%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(m, Phi, c, rho, x, z, u)
   % Linear quadratic problem
   b = -[c - rho * (z - u); -zeros(length(c)-m,1)];
   x = Phi * b;
   x = x(1:length(c));
end