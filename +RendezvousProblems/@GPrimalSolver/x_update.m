%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Inputs:  
% Outputs: - vector x, the update impulsive sequence

% Proximal minimization for Z
function [x] = x_update(indices, p, rho, x, z, u)
    % Impulses update
    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Fuel consumption minimization
        switch (p)
            case 'L1'
                x(sel) = l1_shrinkage(z(sel) - u(sel), 1/rho);
            case 'L2'
                x(sel) = l2_shrinkage(z(sel) - u(sel), 1/rho);
            case 'Linfty'
                x(sel) = lifty_shrinkage(z(sel) - u(sel), 1/rho);
        end
    end
end

%% Auxiliary functions
% Proximal minimization of the L1 norm
function z = l1_shrinkage(x, kappa)
    z = x;
    for i = 1:length(x)
        z(i) = max(0, x(i)-kappa) - max(0, -x(i)-kappa);
    end
end

% Proximal minimization of the L2 norm
function z = l2_shrinkage(x, kappa)
    z = max(0, 1 - kappa/norm(x)) * x;
end

% Proximal minimization of the Lifty norm
function z = lifty_shrinkage(x, kappa)
    z = x - kappa * l1_bproj(x/kappa,1);
end