%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Inputs:  
% Outputs: - vector z, the update impulsive sequence

% Proximal minimization for X
function [x] = x_update(Phi, b, rho, x, z, u)
    % Impulses update
    e = [rho * (z - u); b]; 
    sol = Phi * e;

    x = sol(1:length(x));
end