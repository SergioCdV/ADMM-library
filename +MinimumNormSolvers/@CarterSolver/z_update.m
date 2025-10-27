%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(n, p, q, ~, umax, K, Phi, b, rho, x, ~, u)
    % Pre-allocation 
    y = x + u;

    % Fuel consumption minimization
    switch (p)
        case src.VectorNorm.L1
            hand_ = @(z)src.MinL1Prox.projection( 1/rho, z );

        case src.VectorNorm.L2
            hand_ = @(z)src.MinL2Prox.projection( 1/rho, z );

        case src.VectorNorm.Linfty
            hand_ = @(z)src.MinLinftyProx.projection( 1/rho, z );
    end

    % Maximum control ball projection
    if ( umax ~= Inf )
        switch q
            case src.VectorNorm.L1
                proj_handle_ = @(z)src.L1BallProx.projection( umax, z );

            case src.VectorNorm.L2
                proj_handle_ = @(z)src.L2BallProx.projection( umax, z );

            case src.VectorNorm.Linfty
                proj_handle_ = @(z)src.LinftyBallProx.projection( umax, z );
        end
    else
        proj_handle_ = @(z)( z );
    end

    % Fuel consumption minimization
    idx = 1:n;
    dV = reshape( y(idx), n, []);
    z(:,idx) = hand_( dV );

    % Control authority (this should be parallel projections really)
    z(:,idx) = proj_handle_( z(:,idx) );

    % Cardinality constraint
    cost = p.ComputeVectorNorm( z(:,idx) );

    if ( K ~= Inf )
        [~, pos] = sort( cost, 'descend' );
        sel = pos(K+1:end);
        z(:,sel) = zeros(n, length(sel));
        cost(sel) = zeros(1, length(sel));
    end

    % Primer vector updates
    for i = 1:size(z,2)
        idx = 1 + n * (i - 1) : n * i;
        Phi(:,idx) = -cost * Phi(:,idx);
    end

    z = reshape(z, [], 1);

    pinvA = pinv( Phi );
    Ab = pinvA * b;
    idx = n+1:2*n;
    z(idx) = (eye(n) - pinvA * Phi) * y(idx) + Ab;
end