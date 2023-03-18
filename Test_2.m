%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise II %% 
% Minimize 1/2*sum(huber(A*x - b)) % 

close; 
clear; 
clc

rand('seed', 0);
randn('seed', 0);

% Problem 
m = 5000;       % number of examples
n = 200;        % number of features

x0 = randn(n,1);
A = randn(m,n);
A = A*spdiags(1./normalize(A)',0,n,n); % normalize columns
b = A*x0 + sqrt(0.01)*randn(m,1);
b = b + 10*sprand(m,1,200/m);      % add sparse, large noise

Atb = A'*b;
[L, U] = factor(A);

% Constraint 
A = eye(size(A,2)); 
B = -A;        
c = b;
rho = 1;

% Create the functions to be solved 
objective = @(x)norm(x,1);
X_update = @(x,z,u)(U \ (L \ (Atb + A'*(z - u)) ));
Z_update = @(x,z,u)(rho/(1+rho) * (x+u-b) + (max(0, x+u-b-1-1/rho) - max(0, -x-u+b-1-1/rho))/(1+rho) );

Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);

Problem.alpha = 1;

[x,z,Output] = Problem.solver();