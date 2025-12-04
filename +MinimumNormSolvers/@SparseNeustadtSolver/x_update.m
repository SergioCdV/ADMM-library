%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 04/12/25
% File: x_update.m 
% Issue: 0 
% Validated: 

%% X update %% 
% ADMM problem function to update the X sequence via proximal minimization

function [x] = x_update(Fv, V, c, nx, rho, ~, z, u)
   % Linear problem
   res = z - u;                            % Ax - z residual
   Z = res( 1:nx );                        % Z partition
   Y = res( nx + 1 :end );                 % Y partition
   x = Fv .* ( -c / rho + Z + V * Y );     % Projection onto the Ax - z constraint while minimizing c^Tx
end