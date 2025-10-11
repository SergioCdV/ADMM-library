%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: HalfSpaceProx.m 
% Issue: 0 
% Validated:

%% Half-space Proximal operator %%
% This function defines the proximal operator of a half space

classdef HalfSpaceProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = HalfSpaceProx(a, b, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)HalfSpaceProx.projection(a, b, x), rho )
        end
    end

    methods (Static)
        % Projection onto the halfspace 
        function [y] = projection(a, b, x)
            y = x - max( dot(a,x) - b, 0) / dot(a,a) * a;
        end
    end
end