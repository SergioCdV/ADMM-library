%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the Lasso problem as a proxy to rendezvous

close; 
clear; 
clc

rand('seed', 0);
randn('seed', 0);

n = 100;

x0 = ones(n,1);
for j = 1:3
    idx = randsample(n,1);
    k = randsample(1:10,1);
    x0(ceil(idx/2):idx) = k*x0(ceil(idx/2):idx);
end
b = x0 + randn(n,1);

lambda = 5;

e = ones(n,1);
D = spdiags([e -e], 0:1, n,n);

I = speye(n);
DtD = D'*D;

rho = 1;

% Create the functions to be solved 
objective = @(x,z)(.5*norm(x - b)^2 + lambda*norm(z,1));
X_update = @(x,z,u)((I + rho*DtD) \ (b + rho*D'*(z-u)));
Z_update = @(x,z,u)(max(0, x+u-lambda/rho) - max(0, -x-u-lambda/rho));

% Constraint 
A = D;
B = -eye(n);        
c = zeros(n,1);

% Problem
Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

% Optimization
[x,z,Output] = Problem.solver();

%% Auxiliary functions 
function [x] = x_update(A,R,b,rho,x,z,u)
    x = R \ (R' \ (A'*(b + z - u)));
end