%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: L1BallProx.m 
% Issue: 0 
% Validated:

%% L1-ball Proximal operator %%
% This function defines the proximal operator of an L1-norm ball

classdef L1BallProx < Methods.ProxOperator
    methods
        % Constructor 
        function [obj] = L1BallProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@Methods.ProxOperator( @(x)L1BallProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto an L1-ball 
        function [y] = projection(a, x)
            if ( sum( abs(x) ) > a )
                u = sort(x, 'descend');
                K = 1:length(x); 

                index = cumsum(u) ./ K < u;
                index(1) = 1;

                u = u( logical(index) );
                rho = sum(u - a) / length(u);
                y = max(x - rho, 0);

            else
                y = x;
            end
        end
    end
end