%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise V %% 
% Solve for the minimum of (1/2)*x'*P*x + q'*x + r subject to lb <= x <= ub

close; 
clear; 
clc

randn('state', 0);
rand('state', 0);

n = 100;

% generate a well-conditioned positive definite matrix
% (for faster convergence)
P = rand(n);
P = P + P';
[V, D] = eig(P);
P = V*diag(1+rand(n,1))*V';

q = randn(n,1);
r = randn(1);

l = randn(n,1);
u = randn(n,1);
lb = min(l,u);
ub = max(l,u);

rho = 1;

% Create the functions to be solved 
objective = @(x,z)(0.5*x'*P*x + q'*x + r);
X_update = @(x,z,u)x_update(P,q,rho,x,z,u);
Z_update = @(x,z,u)(z_update(ub,lb,x,z,u));

% Constraint 
A = eye(n); 
B = -A;        
c = zeros(n,1);

% Problem
Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

% Optimization
[x,z,Output] = Problem.solver();

%% Auxiliary functions 
function [x] = x_update(P,q,rho,x,z,u)
    R = chol(P + rho*eye(size(P,1)));
    x = R \ (R' \ (rho*(z - u) - q));
end

function [z] = z_update(ub,lb,x,z,u)
 z = min(ub, max(lb, x + u));
end