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
            y = x;
            norm_ = vecnorm( y, 2, 1 );
            idx = norm_ > a;
            if any(idx)
                y(:,idx) = a * x(:,idx) ./ norm_(idx);
            end
        end
    end
end