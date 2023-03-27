%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the CW optimal L2 problem using ADMM and PVT

close; 
clear; 
clc

%% Define the rendezvous problem and the STM 
% Time span
t = 0:2e-3:pi; 
n = 1;

% Relative initial conditions 
x0 = [1; 2; 3; -1; -2*n*1; 2]; 

% Relative final conditions 
xf = zeros(6,1);

% CW out-of-plane STM
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

% Residual or missvector
b = xf-M*x0;

%% Define the ADMM problem 
rho = 10;        % AL parameter 

% pre-factor
Atb = pinv(Phi)*b;
pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);

p = repmat(3,1,length(t));
cum_part = cumsum(p);

% Create the functions to be solved 
Obj = @(x,z)(objective(cum_part, z));
X_update = @(x,z,u)(x_update(pInvA, Atb, x, z, u));
Z_update = @(x,z,u)(z_update(rho, p, x, z, u));

% Constraint 
[m,n] = size(Phi); 
A = eye(n);
B = -eye(m);        
c = zeros(m,1);

% Problem
Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

%% Optimization
[x, z, Output] = Problem.solver();

%% Apply the pruner 
dV = reshape(x(:,end), 3, []);

cost_admm = sum(sqrt(dot(dV,dV,1)));
[dV, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV);

%% Outcome 
res(:,1) = b-Phi*x(:,end);
res(:,2) = b-Phi*reshape(dV, [], 1);
ratio = 1-cost/cost_admm;

%% Results 
figure
plot(1:Output.Iterations, Output.objval); 
ylabel('Fval')
xlabel('Iterations')
grid on;

figure
plot(t, dV); 
grid on;
ylabel('dV')
xlabel('Time')

figure
stem(t, sqrt(dot(dV,dV,1)), 'filled'); 
grid on;
ylabel('dV')
xlabel('Time')

%% Auxiliary functions 
function [x] = x_update(pInvA, Atb, x, z, u) 
   x = pInvA*(z-u)+Atb;
end

function [z] = z_update(rho, p, x, z, u)
    % cumulative partition
    cum_part = cumsum(p);

    % z-update
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:cum_part(i);
        z(sel) = shrinkage(x(sel) + u(sel), 1/rho);
        start_ind = cum_part(i) + 1;
    end
end

function z = shrinkage(x, kappa)
    z = max(0, x-kappa) - max(0, -x-kappa);
end

function [p] = objective(cum_part, z)
    obj = 0;
    start_ind = 1;
    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj = obj + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end
    p = ( obj );
end