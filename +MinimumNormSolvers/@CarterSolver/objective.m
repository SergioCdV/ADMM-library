%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: objective.m 
% Issue: 0 
% Validated: 

%% Objective %% 
% ADMM problem function to compute the fuel objective

function [cost] = objective(p, z)
    % Initialization 
    cost = sum( p.ComputeVectorNorm(z) );
end