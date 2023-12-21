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

function [p] = objective(p, cum_part, z)
    % Initialization 
    obj = 0;
    start_ind = 1;

    switch (p)
        case 'L1'
            for i = 1:length(cum_part)
                sel = start_ind:cum_part(i);
                obj = obj + sum( abs(z(sel)) );
                start_ind = cum_part(i) + 1;
            end

        case 'L2'
            for i = 1:length(cum_part)
                sel = start_ind:cum_part(i);
                obj = obj + norm(z(sel));
                start_ind = cum_part(i) + 1;
            end

        case 'Linfty'
            for i = 1:length(cum_part)
                sel = start_ind:cum_part(i);
                obj = obj + norm(z(sel), "inf");
                start_ind = cum_part(i) + 1;
            end
    end

    p = obj;
end