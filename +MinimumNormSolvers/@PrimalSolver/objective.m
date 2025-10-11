%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: objective.m 
% Issue: 0 
% Validated: 

%% Objective %% 
% ADMM problem function to compute the fuel objective

function [cost] = objective(p, cum_part, z)
    % Initialization 
    cost = 0;
    start_ind = 1;

    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        cost = cost + p.ComputeVectorNorm( z(sel) );
        start_ind = cum_part(i) + 1;
    end
end