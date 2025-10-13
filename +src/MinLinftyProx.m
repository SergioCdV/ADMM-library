%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: MinLinftyProx.m 
% Issue: 0 
% Validated:

%% Minimum Linfty Proximal operator %%
% This function defines the proximal operator for the minimum Linfty norm

classdef MinLinftyProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = MinLinftyProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)MinLinftyProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto the minimum Linfty norm 
        function [y] = projection(kappa, x)
            y = x - kappa * src.MinL1Prox.projection(1, x / kappa);
        end
    end
end