%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: z_update.m 
% Issue: 0 
% Validated: 

%% Z update %% 
% ADMM problem function to update the Z sequence via proximal minimization

function [z] = z_update(n, og_idx, sigma_map, sigma_unique, q, Phi, b, ~, x, ~, u)
    % Constants
    m = size(Phi,2);
    y = x + u;

    % Projection of the primer vector on the unit+sigma lq-ball
    switch (q)
        case src.VectorNorm.L1
            handl_ = @(p, sigma)src.L1EpigraphProx.projection(1, p, sigma);

        case src.VectorNorm.L2
            handl_ = @(p, sigma)src.LorentzProx.projection(p, sigma);

        case src.VectorNorm.Linfty
            handl_ = @(p, sigma)src.LinfEpigraphProx.projection(1, p, sigma);
    end

    % Primer vectors
    Tk = length(og_idx);               % Number of Lagrange multipliers

    primer_idx = og_idx(end) + Tk + 1 : length(x);

    p = y(primer_idx);                 % Primer vector
    p = reshape(p, n, []);

    % Vectorization of variables
    t = y(Tk + og_idx).';              % Lagrange multipliers associated to the control bound
    t = t(sigma_map);                  % Lagrange multiplier associated to each impulse

    % Projection onto the corresponding epigraph
    [p, t] = handl_( p, t );

    % Lagrange multiplier update (projection on a half space)
    lambda = y(1:m);
    lambda = src.HalfSpaceProx.projection(b, 0, lambda);

    % Projection onto the positive orthant
    sigma = y(og_idx);
    sigma = max(sigma, 0);

    % Final vector
    z = [lambda; sigma; t(sigma_unique).'; reshape(p, [], 1)];
end