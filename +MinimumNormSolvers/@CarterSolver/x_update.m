%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(n, q, pInvA, Atb, x, z, u) 
    N = length(x) / 2;
    y = z(1:N) - u(1:N);
   
    % Impulses update (proximal minimization of the flow indicator function: Ax = b)
    x(1:N) = pInvA * y + Atb;         

    % Primer vector update (projection onto the unit q-ball)
    dV = reshape(y, n, []);
    p = dV;
    
    p_norm   = q.ComputeVectorNorm( dV );
    idx      = p_norm ~= 0;
    p(:,idx) = dV(:,idx) ./ p_norm(idx);

    x(N+1:end) = reshape(p, [], 1);
end