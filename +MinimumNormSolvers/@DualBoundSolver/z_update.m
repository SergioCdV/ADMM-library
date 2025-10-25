%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(indices, tw_idx, q, Phi, b, ~, x, ~, u)
    % Constants 
    m = size(Phi,2);
    y = x + u;

    % Projection of the primer vector on the unit+sigma lq-ball
    switch (q)
        case src.VectorNorm.L1
            handl_ = @(p, sigma)src.L1EpigraphProx.projection(1, p, sigma);

        case src.VectorNorm.L2
            handl_ = @(p, sigma)src.L2EpigraphProx.projection(1, p, sigma);

        case src.VectorNorm.Linfty
            handl_ = @(p, sigma)src.LinfEpigraphProx.projection(1, p, sigma);
    end

    % Primer vector update
    N = max(tw_idx);            % Number of discretised time windows
    p = y(m+N+1:end);           % Primer vector
    sigma = y(m+1:m+N);         % Lagrange multipliers associated to the control bound

    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);
        [p(sel), sigma( tw_idx(i) )] = handl_( p(sel), sigma( tw_idx(i) ) );
        start_ind = indices(i) + 1;
    end

    % Lagrange multiplier update (projection on a half space)
    lambda = y(1:m);
    lambda = src.HalfSpaceProx.projection(b, 0, lambda);

    % Final vector
    z = [lambda; sigma; p];
end