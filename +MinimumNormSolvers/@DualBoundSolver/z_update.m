%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(m, n, Nopp, og_idx, q, b, ~, x, ~, u)
    % Constants
    y = x + u;

    % Lagrange multiplier update (projection on a half space)
    lambda = y(1:m);
    lambda = src.HalfSpaceProx.projection(b, 0, lambda);

    % Projection onto the positive orthant of the Lagrange multipliers associated to the control bound
    Tk = length(og_idx);               % Number of Lagrange multipliers
    sigma = y(og_idx);
    t = y(Tk + og_idx);                
    
    sigma = max(sigma, 0);
    t = max(t, 0);

    % Window bounds 
    tj_idx = og_idx(end) + Tk + 1 : og_idx(end) + Tk + Nopp;
    tj = y(tj_idx);

    if ( size(tj,1) ~= 1 )
        tj = tj.';
    end

    % Primer vectors
    primer_idx = tj_idx(end) + 1 : length(x);

    p = y(primer_idx);                 
    p = reshape(p, n, []);

    % Projection of the primer vector on the t_j lq-ball
    switch (q)
        case src.VectorNorm.L1
            handl_ = @(p, sigma)src.L1EpigraphProx.projection(1, p, sigma);

        case src.VectorNorm.L2
            handl_ = @(p, sigma)src.LorentzProx.projection(p, sigma);

        case src.VectorNorm.Linfty
            handl_ = @(p, sigma)src.LinfEpigraphProx.projection(1, p, sigma);
    end

    % Projection onto the corresponding epigraph
    [p, tj] = handl_( p, tj );

    % Final vector
    z = [lambda; sigma; t; reshape(tj, [], 1); reshape(p, [], 1)];
end