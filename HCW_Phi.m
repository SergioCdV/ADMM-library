%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 24/09/23
% File: HCW_Phi.m 
% Issue: 0 
% Validated: 

%% Fundamental solution matrix of HCW %% 
% Implementation of the HCW STM

function [STM] = HCW_Phi(n, t)
    c = cos(n*t); 
    s = sin(n*t); 

    STM = [1 0 6*(t-s) 4/n*s-3*t 0 2/n*(1-c); ...
           0 c 0 0 s/n 0; ...
           0 0 4-3*c 2/n*(c-1) 0 s/n; ...
           0 0 6*n*(1-c) 4*c-3 0 2*s; ...
           0 -n*s 0 0 c 0; ...
           0 0 3*n*s -2*s 0 c];
end