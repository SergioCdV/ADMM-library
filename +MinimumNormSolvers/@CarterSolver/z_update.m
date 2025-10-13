%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(indices, p, q, umin, umax, K, Phi, b, rho, x, z, u)
    % Pre-allocation 
    y = x + u;

    % Impulses update
    n = length(x) / 2;

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
    if (umax ~= Inf)
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

    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Fuel consumption minimization
        z(sel) = hand_( 1/rho, y(sel) );

        % Maximum control ball projection 
        z(sel) = proj_handle_( z(sel) );

        % Maximum control ball projection   
        norm = q.ComputeVectorNorm( z(sel) ); 
        Phi(:,sel) = -norm * Phi(:,sel);
                   
        start_ind = indices(i) + 1;
    end

    % Cardinality constraint
    if (K ~= Inf)
        dV = reshape(z, indices(1), []);
        cost = p.ComputeVectorNorm( dV );
    
        [~, pos] = sort( cost, 'descend');
        
        index = pos(K+1:end);
        for i = 1:length(index)
            z(1 + indices(1) * (index(i)-1): indices(1) * index(i)) = zeros(indices(1), 1);
        end
    end

    % Primer vector updates
    pinvA = pinv(Phi);
    Ab = pinvA * b;
    idx = n+1:2*n;
    z(idx) = (eye(n) - pinvA * Phi) * y(idx) + Ab;
end