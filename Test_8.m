%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the Lasso problem

close; 
clear; 
clc

randn('seed', 0);
rand('seed', 0);

m = 1500;                           % amount of data
K = 200;                            % number of blocks
partition = randi(50, [K 1]);
n = sum(partition);                 % number of features
p = 100/n;                          % sparsity density

% generate block sparse solution vector
x = zeros(n,1);
start_ind = 1;
cum_part = cumsum(partition);
for i = 1:K
    x(start_ind:cum_part(i)) = 0;
    if( rand() < p)
        % fill nonzeros
        x(start_ind:cum_part(i)) = randn(partition(i),1);
    end
    start_ind = cum_part(i)+1;
end

% generate random data matrix
A = randn(m,n);

% normalize columns of A
A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns

% generate measurement b with noise
b = A*x + sqrt(0.001)*randn(m,1);

rho = 1;

% lambda max
start_ind = 1;
for i = 1:K
    sel = start_ind:cum_part(i);
    lambdas(i) = norm(A(:,sel)'*b);
    start_ind = cum_part(i) + 1;
end
lambda_max = max(lambdas);

% regularization parameter
lambda = 0.1*lambda_max;

xtrue = x;   % save solution

% pre-factor
[L, U] = factor(A,rho);
Atb = A'*b;

% cumulative partition
cum_part = cumsum(p);

% Create the functions to be solved 
Obj = @(x,z)(objective(A, b, lambda, cum_part, x, z));
X_update = @(x,z,u)(x_update(A, Atb, L, U, rho, x, z, u));
Z_update = @(x,z,u)(z_update(lambda, rho, p,x,z,u));

% Constraint 
[m,n] = size(A); 
A = eye(n);
B = -eye(n);        
c = zeros(n,1);

% Problem
Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
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

function [x] = x_update(A, Atb, L, U, rho, x, z, u)
    [m,n] = size(A); 

    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )              % if skinny
       x = U \ (L \ q);
    else                      % if fat
       x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end
end

function [z] = z_update(lambda, rho, p, x, z, u)
    % cumulative partition
    cum_part = cumsum(p);

    % z-update
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:cum_part(i);
        z(sel) = shrinkage(x(sel) + u(sel), lambda/rho);
        start_ind = cum_part(i) + 1;
    end
end

function z = shrinkage(x, kappa)
    z = max(0,1 - kappa/norm(x))*x;
end

function [p] = objective(A, b, lambda, cum_part, x, z)
    obj = 0;
    start_ind = 1;
    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj = obj + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end
    p = ( 1/2*sum((A*x - b).^2) + lambda*obj );
end