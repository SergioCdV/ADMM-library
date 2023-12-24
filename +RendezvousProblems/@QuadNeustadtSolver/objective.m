%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: objective.m 
% Issue: 0 
% Validated: 

%% Objective %% 
% ADMM problem function to compute the fuel objective

% Inputs:  
% Outputs: - vector z, the update impulsive sequence

function [p] = objective(c, x)
    p = -dot(c,x); 
end