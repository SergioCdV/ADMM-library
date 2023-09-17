%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Inputs:  
% Outputs: - vector z, the update impulsive sequence

function [z] = z_update(indices, umin, umax, K, rho, x, z, u) 

    % Fuel consumption L1 minimization
    for i = 1:length(z)
        z(i) = l1_shrinkage(x(i) + u(i), 1/rho);
    end

    start_ind = 1;
    if (umax ~= Inf)
        for i = 1:length(indices)
            sel = start_ind:indices(i);
    
            % Maximum control ball projection
            z(sel) = lifty_bproj(z(sel), umax);

            start_ind = indices(i) + 1;
        end
    end

    % Cardinality constraint
    if (K ~= Inf)
        dV = reshape(z, indices(1), []);
        cost = sum( abs(dV), 1);
    
        [~, pos] = sort( cost, 'descend');
        
        index = pos(K+1:end);
        for i = 1:length(index)
            z(1 + indices(1) * (index(i)-1): indices(1) * index(i)) = zeros(indices(1), 1);
        end
    end
end

%% Auxiliary functions 
% Proximal minimization of the L1 norm
function z = l1_shrinkage(x, kappa)
    z = max(0, x-kappa) - max(0, -x-kappa);
end

% Projection onto the Lifty ball 
function [x] = lifty_bproj(x, a)
    for i = 1:length(x)
        if (x(i) < -a)
            x(i) = -a;
        elseif (x(i) > a)
            x(i) = a;
        end
    end
end