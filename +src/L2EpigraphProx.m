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
            ep = a + t;

            % Inactive constraint
%             zero_idx = sigma <= 0;
%             y(:,zero_idx) = src.L2BallProx.projection(a, x(:,zero_idx));
            
            % Active constraints
            norm_ = vecnorm( x, 2, 1 );
            idx = norm_ > ep;
            if ( sum(idx) )
                r = norm_(idx);
                C = r + ep(idx);
                alpha = C ./ (2 * r);
                alpha = alpha .* sign(alpha);
                y(:,idx) = alpha .* x(:,idx);
                sigma(idx) = alpha .* r - a;                         
            end
        end
    end
end