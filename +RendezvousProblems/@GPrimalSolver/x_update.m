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

function [x] = x_update(indices, p, q, umin, umax, K, pInvA, Atb, rho, x, z, u) 
   x = pInvA * (z-u) + Atb;                     % Impulses update (proximal minimization of the flow indicator function: Ax = b)

   % Maximum control ball projection
    if (umax ~= Inf)
        start_ind = 1;
        for i = 1:length(indices)
            sel = start_ind:indices(i);
            switch (q)
                case 'L1'
                    x(sel) = l1_bproj(x(sel), umax);
                case 'L2'
                    x(sel) = l2_bproj(x(sel), umax);
                case 'Linfty'
                    x(sel) = lifty_bproj(x(sel), umax);
            end
            start_ind = indices(i) + 1;
        end
    end

    % Cardinality constraint
    if (K ~= Inf)
        dV = reshape(x, indices(1), []);
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
            x(1 + indices(1) * (index(i)-1): indices(1) * index(i)) = zeros(indices(1), 1);
        end
    end
end

%% Auxiliary functions
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