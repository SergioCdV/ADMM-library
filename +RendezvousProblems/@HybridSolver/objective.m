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
%          - vector x, the update impulsive sequence

function [p] = objective(cum_part, x, z)
    % Initialization 
    obj1 = 0;
    obj2 = 0;
    start_ind = 1;

    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj1 = obj1 + sum( abs(x(sel)) );
        start_ind = cum_part(i) + 1;
    end

    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj1 = obj1 + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end

    p = obj1 + obj2;
end