%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(n, ~, umax, K, rho, x, ~, u) 
    % Pre-allocation 
    y = x + u;
    y = reshape(y, n, []); 

    % Fuel consumption L1 minimization
    z = src.MinL1Prox.projection( 1/rho, y );

    % Maximum control ball projection
    z = Methods.LinfBallProx.projection( umax, z );

    % Cardinality constraint
    if ( K ~= Inf )
        [~, pos] = sort( cost, 'descend' );
        sel = pos(K+1:end);
        z(:,sel) = zeros(n, length(sel));
    end

    z = reshape(z, [], 1);
end