%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 02/09/23
% File: YA_Phi.m 
% Issue: 0 
% Validated: 

%% Fundamental solution matrix of Y-A %% 
% Implementation of the YA fundamental matrix 

function [STM] = YA_invPhi(mu, h, e, dt, f)
    % Define Carter's constants 
    k = 1 + e * cos(f);

    % Assemble the matrix 
    c = k * cos(f); 
    s = k * sin(f);

    STM = [1-e^2 0 3*e*s*(1/k+1/k^2) -e*s*(1+1/k) 0 -e*c+2; ...
           0 0 0 0 0 0; ...
           0 0 -3*s*(1/k+e^2/k^2) s*(1+1/k) 0 c-2*e; ...
           0 0 -3*(c/k+e) c*(1+1/k) 0 -s; ...
           0 0 0 0 0 0; ...
           0 0 3*k+e^2-1 -k^2 0 e*s]/(1-e^2);
end