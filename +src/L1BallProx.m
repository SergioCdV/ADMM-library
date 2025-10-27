%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: L1BallProx.m 
% Issue: 0 
% Validated:

%% L1-ball Proximal operator %%
% This function defines the proximal operator of an L1-norm ball

classdef L1BallProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = L1BallProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)L1BallProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto an L1-ball 
        function [x] = projection(a, x)
            norm_ = vecnorm( x, 1, 1 );
            idx = norm_ > a;

            if ( any(idx) )
                N = sum(idx);
                m = size(x,1);
                cols = 1:m;
                K = repmat(cols.', 1, N);

                u = sort( x(:,idx), 'descend' );
                U = cumsum(u, 1) - a;                    
                res = U ./ K;

                index = res < u;
                [~, idx_from_bottom] = max( flipud(index), [], 1 );
                last_k = m - idx_from_bottom + 1;
                rho = res(last_k);             

                x(:,idx) = max( x(:,idx) - rho, 0 );
            end
        end
    end
end