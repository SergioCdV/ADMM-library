%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(m, n, nx, q, c, A, U, b, rho, x, ~, u)
    % Partition of data 
    idx = 1 : nx; 
    Z = x(idx) + u(idx);        % Z partition 

    idx = nx + 1 : length(x);
    Y = x(idx) + u(idx);        % Y partition

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
    p = reshape( Z(m+1:end), n, [] );
    p = handl_(p);

    % Lagrange multiplier update (projection on a half space)
    lambda = Z(1:m);
    lambda = src.HalfSpaceProx.projection(c, 0, lambda);

    % Final partition
    Z = [lambda; reshape(p, [], 1)];

    % Projection onto the sparse linear system between lambda and the primer vector
    Y = Y - U.' * ( ( U * Y - b ) ./ A );

    % Final vector
    z = [Z; Y];
end