%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise I %% 
% Solve for the minimum of ||x||_1 s.t. Ax = b % 

close; 
clear; 
clc

rand('seed', 0);
randn('seed', 0);

% Problem 
n = 30;
m = 10;
A = randn(m,n);

x = sprandn(n, 1, 0.1*n);
b = A*x;

xtrue = x;

AAt = A*A';
P = eye(size(A,2)) - A' * (AAt \ A);
q = A' * (AAt \ b);

% Constraint 
A = eye(size(A,2)); 
B = -A;        
c = zeros(size(A,2),1);
rho = 1;

% Create the functions to be solved 
objective = @(x)norm(x,1);
X_update = @(x,z,u)(P*(z - u) + q);
Z_update = @(x,z,u)(max(0, x+u-1/rho) - max(0, -x-u-1/rho));

Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);

Problem.alpha = 1;

[x,z,Output] = Problem.solver();