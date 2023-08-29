%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Inputs:  
% Outputs: - vector x, the update impulsive sequence

function [x] = x_update(m, Phi, c, rho, x, z, u)
   % Impulses update (proximal minimization of the flow indicator function: Ax = b)
   b = [c + rho * (z - u); zeros(length(c)-m,1)];
   x = Phi * b;
   x = x(1:length(c));
end