%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise IX %% 
% Solve for the CW optimal L2 problem using ADMM and PVT for free time

close; 
clear; 
clc

%% Define the rendezvous problem %%
% Time span
t = linspace(0,1,200);        % Nondimensional
n = 1;                         % Target mean motion

% Final time restrictions
LB = pi;                        % Minimum TOF
UB = 100;                       % Maximum TOF

% Maximum DV 
dVmax = 0.5;                     % Maximum L2 norm of each impulse

% Relative initial conditions 
x0 = [1; 2; 3; -1; -2*n*1; 2]; 

% Relative final conditions 
xf = zeros(6,1);

% Number of impulses
p = repmat(3,1,length(t));
cum_part = cumsum(p);

%% Upper level solver 
% AL parameter 
rho = 1;        % AL parameter 

% Update functions 
Obj = @(x,z)(objective(cum_part, z));
X_update = @(x,z,u)(x_update(t, n, x0, xf, x, z, u));
Z_update = @(x,z,u)(z_update(rho, dVmax, LB, UB, p, x, z, u));

% Dual-primal constraint 
d = 3*length(t)+1; 
A = eye(d);
B = -eye(d);        
c = zeros(d,1);

% Problem
UpperProblem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
UpperProblem.alpha = 1;

%% Optimization 
[x, z, Output] = UpperProblem.solver();

% Solve for the optimal sequence 
dV = reshape(x(1:end-1,end), 3, []);
TOF = x(end);

% Analysis 
cost_admm = sum(sqrt(dot(dV,dV,1)));

% Generate the STM 
t = TOF * t;
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
res(:,1) = b-Phi*x(1:end-1,end);

%% Pruning
% [dV, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV);
% res(:,2) = b-Phi*reshape(dV, [], 1);
% ratio = 1-cost/cost_admm;

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

%% Lower problem update functions
function [x] = x_update(t, n, x0, xf, x, z, u) 
    % Dimensionalization
    t = z(end) * t;

    % Generate the STM 
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

    % Pre-factoring
    Atb = pinv(Phi)*b;
    pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);
    x(1:end-1) = pInvA*(z(1:end-1)-u(1:end-1))+Atb;
    x(end) = z(end)-u(end);
end

function [z] = z_update(rho, dVmax, LB, UB, p, x, z, u)
    % cumulative partition
    cum_part = cumsum(p);

    % Impulses updates
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:cum_part(i);
        z(sel) = shrinkage(x(sel) + u(sel), 1/rho);

        if (norm(z(sel)) > dVmax)
            z(sel) = dVmax * z(sel)/norm(z(sel));
        end

        start_ind = cum_part(i) + 1;
    end

    % TOF update (box projection)
    if (z(end) < LB)
        z(end) = LB;
    elseif (z(end) > UB)
        z(end) = UB;
    end
end

function [z] = shrinkage(x, kappa)
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
