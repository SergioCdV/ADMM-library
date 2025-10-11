%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise IV %% 
% Solve for the minimum of c.'x s.t. Ax = b, x >= 0 % 

close; 
clear; 
clc

randn('state', 0);
rand('state', 0);

n = 500;  % dimension of x
m = 400;  % number of equality constraints

c  = rand(n,1) + 0.5;    % create nonnegative price vector with mean 1
x0 = abs(randn(n,1));    % create random solution vector

A = abs(randn(m,n));     % create random, nonnegative matrix A
b = A*x0;

rho = 1;

% Create the functions to be solved 
objective = @(x,z)(c'*x);
X_update = @(x,z,u)x_update(A,b,c,rho,x,z,u);
Z_update = @(x,z,u)(z_update(x,z,u));

% Constraint 
A = eye(size(A,2)); 
B = -A;        
c = zeros(size(A,2),1);

% Problem
Problem = ADMM_solver(objective, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

% Optimization
[x,z,Output] = Problem.solver();

%% Auxiliary functions 
function [x] = x_update(A,b,c,rho,x,z,u)
    [m,n] = size(A);
    tmp = [ rho*eye(n), A'; A, zeros(m) ] \ [ rho*(z - u) - c; b ];
    x = tmp(1:n,1);
end

function [z] = z_update(x,z,u)
    z = x+u;
    z = max(0,z);
end