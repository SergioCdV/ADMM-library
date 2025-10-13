%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(indices, umin, umax, K, rho, x, z, u) 
    % Fuel consumption L1 minimization
    y = x + u;
    z = src.MinL1Prox.projection( 1/rho, y );

    % Maximum control ball projection
    if (umax ~= Inf)
        start_ind = 1;
        for i = 1:length(indices)
            sel = start_ind:indices(i);
            z(sel) = Methods.LinfBallProx.projection( umax, z(sel) );
            start_ind = indices(i) + 1;
        end
    end

    % Cardinality constraint
    if (K ~= Inf)
        dV = reshape(z, indices(1), []);
        cost = src.VectorNorm.L1.ComputeVectorNorm( dV );
    
        [~, pos] = sort( cost, 'descend');
        
        index = pos(K+1:end);
        for i = 1:length(index)
            z(1 + indices(1) * (index(i)-1): indices(1) * index(i)) = zeros(indices(1), 1);
        end
    end
end