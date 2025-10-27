%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 26/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(n, p, q, ~, umax, K, rho, x, ~, u)
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

    % Impulses update
    y = reshape(y, n, []);

    % Fuel minimization 
    z = hand_( y ); 

    % Control authority (this should be parallel projections really)
    z = proj_handle_( z );

    % Cardinality constraint
    if ( K ~= Inf )
        cost = p.ComputeVectorNorm( z );
        [~, pos] = sort( cost, 'descend' );
        idx = pos(K+1:end);
        z(:,idx) = zeros(n, length(idx));
    end

    z = reshape(z, [], 1);
end