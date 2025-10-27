%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(n, q, Phi, b, ~, x, ~, u)
    % Constants 
    m = size(Phi,2);
    y = x + u;

    % Projection of the primer vector on the unit lq-ball
    switch (q)
        case src.VectorNorm.L1
            handl_ = @(p)src.L1BallProx.projection(1, p);

        case src.VectorNorm.L2
            handl_ = @(p)src.L2BallProx.projection(1, p);

        case src.VectorNorm.Linfty
            handl_ = @(p)src.LinfBallProx.projection(1, p);
    end

    % Primer vector update
    p = reshape( y(m+1:end), n, [] );
    p = handl_(p);

    % Lagrange multiplier update (projection on a half space)
    lambda = y(1:m);
    lambda = src.HalfSpaceProx.projection(b, 0, lambda);

    % Final vector
    z = [lambda; reshape(p, [], 1)];
end