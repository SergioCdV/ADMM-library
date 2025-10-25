%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: L2EpigraphProx.m 
% Issue: 0 
% Validated:

%% L2-epigraph Proximal operator %%
% This function defines the proximal operator of an epigraph of an L2-norm ball

classdef L2EpigraphProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = L2EpigraphProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)L2EpigraphProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto an L2-ball epigraph
        function [y, sigma] = projection(a, x, t)
            % Projection over the epigraph
            y = x; 
            sigma = t; 

            norm_ = vecnorm( x );
            idx = norm_ > a + t;

            if ( sum(idx) )
                alpha = (norm_(idx) + a + t) ./ (2 * norm_(idx));
                y(:,idx) = alpha .* x(:,idx);
                sigma = alpha * norm_(idx) - a;                         % << this is not fully vectorized
                
                % Check the nonnegativity of t
                if ( sigma < 0 )
                    sigma = 0;
                    y(:,idx) = a * x(:,idx) / norm_(idx);
                end
            end
        end
    end
end