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
function [z] = z_update(indices, q, b, rho, x, z, u)
    % Constants 
    m = length(b);

    % Primer vector update
    p = x(m+1:end) + u(m+1:end);

    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Projection of the primer vector on the unit lq-ball
        switch (q)
            case 'L1'
                p(sel) = l1_bproj(p(sel), 1);
            case 'L2'
                p(sel) = l2_bproj(p(sel), 1);
            case 'Linfty'
                p(sel) = lifty_bproj(p(sel), 1);
        end
                   
        start_ind = indices(i) + 1;
    end

    % Lagrange multiplier update (projection on a half space)
    lambda = x(1:m) + u(1:m);
    lambda = half_space(lambda, b, 0);

    % Final vector
    z = [lambda; p];
    z(m+1:end) = p;
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

% Projection onto the halfspace 
function [x] = half_space(x, a, b)
    x = x - max(dot(a,x)-b, 0)/dot(a,a) * a;
end