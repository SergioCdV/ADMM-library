%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(indices, q, Phi, b, rho, x, z, u)
    % Constants 
    m = size(Phi,2);

    % Primer vector update
    p = x(m+1:end) + u(m+1:end);

    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Projection of the primer vector on the unit lq-ball
        switch (q)
            case VectorNorm.L1
                p(sel) = Methods.L1BallProx.projection(1, p(sel));

            case VectorNorm.L2
                p(sel) = Methods.L2BallProx.projection(1, p(sel));

            case VectorNorm.Linf
                p(sel) = Methods.LinfBallProx.projection(1, p(sel));
        end
                   
        start_ind = indices(i) + 1;
    end

    % Lagrange multiplier update (projection on a half space)
    lambda = x(1:m) + u(1:m);
    lambda = Methods.HalfSpaceProx.projection(b, 0, lambda);

    % Final vector
    z = [lambda; p];
    z(m+1:end) = p;
end