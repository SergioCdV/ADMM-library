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

%% Define the rendezvous problem and the STM 
% Time span
t = 0:1e-1:pi; 

% Relative initial conditions 
x0 = [1; 2]; 

% Relative final conditions 
xf = zeros(2,1);

% CW out-of-plane STM
STM = [cos(t.') sin(t.') -sin(t.') cos(t.')]; 
Phi = zeros(2, 2 * length(t));
M = reshape(STM(end,:), [2 2]);
for i = 1:length(t)
    Phi(:,1+size(x0,1)*(i-1):size(x0,1)*i) = M * reshape(STM(i,:), [2 2])^(-1);
end

% Residual or missvector
b = xf-M*x0;

%% Define the ADMM problem 
rho = 1;        % AL parameter 
lambda = 1e5;     % Relative weight of the minimum fuel objective

% pre-factor
[L, U] = factor(Phi,rho);
Atb = Phi'*b;

p = repmat(2,1,length(t));
cum_part = cumsum(p);

% Create the functions to be solved 
Obj = @(x,z)(objective(Phi, b, lambda, cum_part, x, z));
X_update = @(x,z,u)(x_update(Phi, Atb, L, U, rho, x, z, u));
Z_update = @(x,z,u)(z_update(lambda, rho, p,x,z,u));

% Constraint 
[m,n] = size(Phi); 
A = eye(n);
B = -eye(n);        
c = zeros(n,1);

% Problem
Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

%% Optimization
[x, z, Output] = Problem.solver();

%% Outcome 
res = b-Phi*x(:,end);

%% Results 
plot(Output.objval); 
ylabel('Fval')
xlabel('Iterations')

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