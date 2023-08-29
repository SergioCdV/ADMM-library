%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

% Inputs:  
% Outputs: - vector x, the update impulsive sequence

function [x] = x_update(n, q, pInvA, Atb, x, z, u) 
    N = length(x)/2;
   
    % Impulses update (proximal minimization of the flow indicator function: Ax = b)
    x(1:N) = pInvA * (z(1:N)-u(1:N)) + Atb;         

    % Primer vector update (projection onto the unit q-ball)
    dV = reshape(z(1:N)-u(1:N), n, N/n);
    switch (q)
        case 'L1'
            p_norm = sum( abs(dV), 1 );
        case 'L2'
            p_norm = sqrt( dot(dV, dV, 1) );
        case 'Linfty'
            p_norm = max( abs(dV) );
    end

    p = dV;
    p(:,p_norm ~= 0) = dV(:,(p_norm ~= 0)) ./ p_norm(p_norm ~= 0);
    x(N+1:end) = reshape(p, [], 1);
end