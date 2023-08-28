%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

% Inputs:  
% Outputs: - vector z, the update impulsive sequence

% Proximal minimization for Z
function [z] = z_update(indices, p, q, umin, umax, K, rho, x, z, u)
    % Impulses update
    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Fuel consumption minimization
        switch (p)
            case 'L1'
                z(sel) = l1_shrinkage(x(sel) + u(sel), 1/rho);
            case 'L2'
                z(sel) = l2_shrinkage(x(sel) + u(sel), 1/rho);
            case 'Linfty'
                z(sel) = lifty_shrinkage(x(sel) + u(sel), 1/rho);
        end

        % Maximum control ball projection
        if (umax ~= Inf)
            switch (q)
                case 'L1'
                    z(sel) = l1_bproj(z(sel), umax);
                case 'L2'
                    z(sel) = l2_bproj(z(sel), umax);
                case 'Linfty'
                    z(sel) = lifty_bproj(z(sel), umax);
            end
        end
                   
        start_ind = indices(i) + 1;
    end

    % Cardinality constraint
    if (K ~= Inf)
        dV = reshape(z, indices(1), []);
        switch (p)
            case 'L1'
                cost = sum( abs(dV), 1);
            case 'L2'
                cost = sqrt( dot(dV, dV, 1) );
            case 'Linfty'
                cost = max(abs(dV), [], 1);
        end
    
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

% Proximal minimization of the L2 norm
function z = l2_shrinkage(x, kappa)
    z = max(0, 1 - kappa/norm(x)) * x;
end

% Proximal minimization of the Lifty norm
function z = lifty_shrinkage(x, kappa)
    z = max(0, x-kappa) - max(0, -x-kappa);
end

% Projection onto the L1 ball 
function [x] = l1_bproj(x, a)
    if (sum(abs(x)) > a)
        u = sort(x,'descend');
        K = 1:length(x); 

        index = false * ones(1,length(K));
        for i = 1:length(K)
            index(i) = sum(u(1:K(i)))/K(i) < u(K(i));
        end
        index(1) = 1;
        u = u(logical(index));
        rho = sum(u-a)/length(u);
        x = max(x-rho,0);
    end
end

% Projection onto the L2 ball 
function [x] = l2_bproj(x, a)
    Norm = norm(x);
    if (Norm > a)
       x = a * x / Norm;
    end
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