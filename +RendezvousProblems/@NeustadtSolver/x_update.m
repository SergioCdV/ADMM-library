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

function [x] = x_update(m, Phi, c, Phi2, rho, x, z, u)
   % Linear quadratic problem
   b = -[c - rho * (z - u); -zeros(length(c)-m,1)];
   x = Phi * b;
   x = x(1:length(c));

    % Maximization of lambda 
%     lambda = c(1:m)/rho + (z(1:m) - u(1:m));
% 
%     % Solve for the primer vectors
%     p = (eye(size(Phi2,2))-pinv(Phi2)*Phi2) * (z(m+1:end)-u(m+1:end)) + Phi2\lambda;
% 
%     e = x - [lambda;p];
end