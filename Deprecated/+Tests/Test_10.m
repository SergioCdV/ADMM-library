%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise X %% 
% Solve for the parallel projections

close; 
clear; 
clc

% Matrix creation 
% Time span
t = linspace(0,pi,10);
n = 1;

% Relative initial conditions 
x0 = [1; 2; 3; -1; -2*n*1; 2]; 

% Relative final conditions 
xf = zeros(6,1);

% CW STM
STM(:,1:6) = [4-3*cos(n*t.') -6*n*t.'+6*sin(n*t.') zeros(length(t),1) 3*n*sin(n*t.') -6*n+6*n*cos(n*t.') zeros(length(t),1)]; 
STM(:,7:12) = [zeros(length(t),1) ones(length(t),1) zeros(length(t),4)];
STM(:,13:18) = [zeros(length(t),2) cos(n*t.') zeros(length(t),2) -n*sin(n*t.')];
STM(:,19:24) = [sin(n*t.')/n -2/n+2*cos(n*t.')/n zeros(length(t),1) cos(n*t.') -2*sin(n*t.') zeros(length(t),1)];
STM(:,25:30) = [2/n-2*cos(n*t.')/n 4/n*sin(n*t.')-3*t.' zeros(length(t),1) 2*sin(n*t.') -3+4*cos(n*t.') zeros(length(t),1)];
STM(:,31:36) = [zeros(length(t),2) sin(n*t.')/n zeros(length(t),2) cos(n*t.')];

Phi = zeros(6, 3 * length(t));
M = reshape(STM(end,:), [6 6]);
for i = 1:length(t)
    Phi(:,1+3*(i-1):3*i) = M * reshape(STM(i,:), [6 6])^(-1) * [zeros(3); eye(3)];
end

% Create the functions to be solved 
X_update = @(x,z,u)(x_update(Phi, length(t), x, z, u));

% Problem
Problem = parallel_projections(X_update, length(t), 6);

% Optimization
[x,z,Output] = Problem.solver();

%% Auxiliary functions 
function [x] = x_update(Phi, N, x, z, u)
    % Mean X 
    n = size(x,1)/N;
    m = size(Phi,2)/N;
    xm = zeros(n,1);
    for i = 1:N
        xm = xm + x(1+n*(i-1):n*i,1);
    end
    xm = xm/N;

    % Projection onto the unit Euclidean ball
    for i = 1:N 
        res = xm-u(1+n*(i-1):n*i,1); 
        A = Phi(:,1+m*(i-1):m*i).';
        x(1+n*(i-1):n*i) = projection(A,res);
    end
end

function [x] = projection(A, x)
    % Project onto the matrix subspace 
    y = A*x; 

    % Renormalize
    Norm = norm(y);
    if (Norm)
        x = x+rand(size(x));
    end
end