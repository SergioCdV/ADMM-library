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
m = 1500;       % number of examples
n = 5000;       % number of features
p = 100/n;      % sparsity density

x0 = sprandn(n,1,p);
A = randn(m,n);
A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
b = A*x0 + sqrt(0.001)*randn(m,1);

lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max;

rho = 1;

[L, U] = factor(A, rho);
Atb = A'*b;

% Create the functions to be solved 
objective = @(x,z)(( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) ));
X_update = @(x,z,u)x_update(A,b,Atb,L,U,rho,x,z,u);
Z_update = @(x,z,u)(max(0, x+u-lambda/rho) - max(0, -x-u-lambda/rho));

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
function [L, U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

function [x] = x_update(A,b,Atb,L,U,rho,x,z,u)
    [m, n] = size(A);

    q = Atb + rho*(z - u);
    if ( m >= n ) 
        x = U \ (L \ q); 
    else 
        x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2; 
    end
end