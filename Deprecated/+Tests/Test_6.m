%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VI %% 
% Solve for the minimum of ||Ax - b||_1

close; 
clear; 
clc

rand('seed', 0);
randn('seed', 0);

m = 1000;       % number of examples
n = 100;        % number of features

A = randn(m,n);
x0 = 10*randn(n,1);
b = A*x0;
idx = randsample(m,ceil(m/50));
b(idx) = b(idx) + 1e2*randn(size(idx));
R = chol(A'*A);

rho = 1;

% Create the functions to be solved 
objective = @(x,z)(norm(z,1));
X_update = @(x,z,u)x_update(A,R,b,rho,x,z,u);
Z_update = @(x,z,u)(max(0, x+u-b-1/rho) - max(0, -x-u+b-1/rho));

% Constraint 
B = -eye(m);        
c = b;

% Problem
Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

% Optimization
[x,z,Output] = Problem.solver();

%% Auxiliary functions 
function [x] = x_update(A,R,b,rho,x,z,u)
    x = R \ (R' \ (A'*(b + z - u)));
end