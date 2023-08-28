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

function [x] = x_update(pInvA, Atb, x, z, u) 
   x = pInvA * (z-u) + Atb;                     % Impulses update (proximal minimization of the flow indicator function: Ax = b)
end