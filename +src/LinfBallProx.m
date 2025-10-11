%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: LinfBallProx.m 
% Issue: 0 
% Validated:

%% Linf-ball Proximal operator %%
% This function defines the proximal operator of an Linf-norm ball

classdef LinfBallProx < Methods.ProxOperator
    methods
        % Constructor 
        function [obj] = LinfBallProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@Methods.ProxOperator( @(x)LinfBallProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto an Linf-ball 
        function [y] = projection(a, x)
            y = x;
            idx = x < -a; 
            y(idx) = -a * ones(1,sum(idx)); 

            idx = x > a; 
            y(idx) = a * ones(1,sum(idx)); 
        end
    end
end