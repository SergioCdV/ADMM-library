%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: ProxOperator.m 
% Issue: 0 
% Validated: 

%% Proximal operator %%
% This function contains the abstract definition of a proximal operator

classdef (Abstract) ProxOperator
    % Main properties
    properties 
        rho = 1;                % Augmented Lagrangian parameter 
        projection_handle
    end

    methods
        % Constructor 
        function [obj] = ProxOperator(projection, rho)    
            if exists(rho, 'var')
                obj.rho = rho;
            end

            if isa(projection, 'function_handle')
                obj.projection_handle = projection;
            else
                error('A prox. operator requires of a valid function handle as argument');
            end
        end

        % Perform the projection 
        function [x] = project(obj, y)
            x = obj.projection_handle( y, obj.rho  );
        end
    end
end