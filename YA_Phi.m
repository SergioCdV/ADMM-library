%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 02/09/23
% File: YA_Phi.m 
% Issue: 0 
% Validated: 

%% Fundamental solution matrix of Y-A %% 
% Implementation of the YA fundamental matrix 

function [STM] = YA_Phi(mu, h, e, dt, f)
    % Define Carter's constants 
    k = 1 + e * cos(f);
    I = mu^2/h^3*dt;

    % Assemble the matrix 
    c = k * cos(f); 
    s = k * sin(f);
    s_p = cos(f) + e*cos(2*f);
    c_p = -(sin(f)+e*sin(2*f));

    STM = [1 0 -c*(1+1/k) s*(1+1/k) 0 3*k^2*I; ...
           0 c/k 0 0 s/k 0; ...
           0 0 s c 0 2-3*e*s*I; ...
           0 0 2*s 2*c-e 0 3*(1-2*e*s*I); ...
           0 -s/k 0 0 c/k 0; ...
           0 0 s_p c_p 0 -3*e*(s_p*I+s/k^2)];
end