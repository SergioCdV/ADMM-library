%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: MinL1Prox.m 
% Issue: 0 
% Validated:

%% Minimum L1 Proximal operator %%
% This function defines the proximal operator for the minimum L1 norm

classdef MinL1Prox < src.ProxOperator
    methods
        % Constructor 
        function [obj] = MinL1Prox(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)MinL1Prox.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto the minimum L1 norm 
        function [y] = projection(kappa, x)
            y = max(0, x - kappa) - max(0, -x - kappa);
        end
    end
end