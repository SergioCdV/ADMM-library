%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: MinL2Prox.m 
% Issue: 0 
% Validated:

%% Minimum L2 Proximal operator %%
% This function defines the proximal operator for the minimum L2 norm

classdef MinL2Prox < src.ProxOperator
    methods
        % Constructor 
        function [obj] = MinL2Prox(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)MinL2Prox.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto the minimum L2 norm 
        function [y] = projection(kappa, x)
            norm_x = vecnorm( x );
            idx = norm_x ~= 0;
            y = x; 
            y(:,idx) = max( 0, 1 - kappa ./ norm_x(idx) ) .* x(:,idx);
        end
    end
end