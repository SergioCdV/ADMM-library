%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: L2BallProx.m 
% Issue: 0 
% Validated:

%% L2-ball Proximal operator %%
% This function defines the proximal operator of an L2-norm ball

classdef L2BallProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = L2BallProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)L2BallProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto an L2-ball 
        function [y] = projection(a, x)
            norm = sqrt( dot(x, x, 1) );
            y = x;
            idx = norm ~= 0 & norm > a;
            
            if ( ~isempty(norm(idx)) )
                u = x(:,idx) ./ norm(idx);
                y(:,idx) =  a * u;
            end
        end
    end
end