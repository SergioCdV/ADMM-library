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
t = 0:1e-2:pi; 

% Relative initial conditions 
x0 = [1; 2]; 

% Relative final conditions 
xf = zeros(2,1);

% CW out-of-plane STM
STM = [cos(t.') sin(t.') -sin(t.') cos(t.')]; 
Phi = zeros(2, 1 * length(t));
M = reshape(STM(end,:), [2 2]);
for i = 1:length(t)
    Phi(:,i) = M * reshape(STM(i,:), [2 2])^(-1) * [0; 1];
end

% Residual or missvector
b = xf-M*x0;

%% Define the ADMM problem 
rho = 1;        % AL parameter 

% pre-factor
Atb = pinv(Phi)*b;
pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);

p = repmat(1,1,length(t));
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
dV = reshape(x(:,end), 1, []);

cost_admm = sum(sqrt(dot(dV,dV,1)));
[dV, cost] = PVT_pruner(STM, [0;1], dV);

%% Outcome 
res = b-Phi*x(:,end);
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
plot(t, sqrt(dot(dV,dV,1))); 
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
    z = max(0, 1 - kappa/norm(x))*x;
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