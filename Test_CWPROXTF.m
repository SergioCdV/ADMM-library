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
t = linspace(0,1,20);        % Nondimensional
n = 1;                         % Target mean motion

% Final time restrictions
LB = 0;                        % Minimum TOF
UB = pi;                       % Maximum TOF

% Maximum DV 
dVmax = 3;                     % Maximum L2 norm of each impulse

% Relative initial conditions 
x0 = [1; 2; 3; -1; -2*n*1; 2]; 

% Relative final conditions 
xf = zeros(6,1);

%% Dimensionalization 
LB = LB/UB; 
n = n/UB;
dVmax = dVmax * UB;
x0(4:6) = x0(4:6) * UB;
xf(4:6) = xf(4:6) * LB;
UB = UB/UB;

%% Upper level solver 
% AL parameter 
rho = 1;        % AL parameter 

% Update functions 
Obj = @(x,z)(upper_objective(LB, UB, x, z));
X_update = @(x,z,u)(upperx_update(t, n, x0, xf, dVmax, x, z, u));
Z_update = @(x,z,u)(upperz_update(LB, UB, x, z, u));

% Dual-primal constraint 
n = length(t); 
A = eye(n);
B = -eye(n,1);        
c = zeros(n,1);

% Problem
UpperProblem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
UpperProblem.alpha = 1;

%% Optimization 
[x, z, Output] = UpperProblem.solver();

%% Final sequence solution
 % Generate the STM 
t = z * t;
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

p = repmat(3,1,length(t));
cum_part = cumsum(p);

% Define the problem 
rho = 1;        % AL parameter 

% Update functions 
Obj = @(x,z)(lower_objective(cum_part, z));
X_update = @(x,z,u)(lowerx_update(pInvA, Atb, x, z, u));
Z_update = @(x,z,u)(lowerz_update(rho, p, x, z, u));

% Dual-primal constraint 
[m,n] = size(Phi); 
A = eye(n);
B = -eye(n);        
c = zeros(m,1);

% Problem
LowerProblem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
LowerProblem.alpha = 1;
LowerProblem.QUIET = false;

% Solve for the optimal sequence 
[dV, ~, cost_admm] = LowerProblem.solver();
dV = reshape(dV(:,end), 3, []);
[dV, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV);

%% Analysis 
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

%% Lower problem update functions
function [x] = lowerx_update(pInvA, Atb, x, z, u) 
   x = pInvA*(z-u)+Atb;
end

function [z] = lowerz_update(rho, p, x, z, u)
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

function [p] = lower_objective(cum_part, z)
    obj = 0;
    start_ind = 1;
    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj = obj + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end
    p = ( obj );
end

%% Upper problem functions
function [x] = upperx_update(t, n, x0, xf, dVmax, x, z, u)
    % Dimensionalization
    n = z * n;
    dVmax = dVmax / z;
    x0(4:6) = x0(4:6) / z;
    xf(4:6) = xf(4:6) / z;

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

    p = repmat(3,1,length(t));
    cum_part = cumsum(p);

    % Define the problem 
    rho = 1;        % AL parameter 

    % Update functions 
    Obj = @(x,z)(lower_objective(cum_part, z));
    X_update = @(x,z,u)(lowerx_update(pInvA, Atb, x, z, u));
    Z_update = @(x,z,u)(lowerz_update(rho, p, x, z, u));
    
    % Dual-primal constraint 
    [m,n] = size(Phi); 
    A = eye(n);
    B = -eye(n);        
    c = zeros(m,1);
    
    % Problem
    LowerProblem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
    LowerProblem.alpha = 1;
    LowerProblem.QUIET = false;

    % Solve for the optimal sequence 
    [dV, ~, ~] = LowerProblem.solver();
    dV = reshape(dV(:,end), 3, []);

    % Individual projection over each set (Eucliden ball)
    for i = 1:size(dV,2)
        index = norm(dV(:,i)) > dVmax;
        x(i,1) = dVmax * x(i,1) * (index) + x(i,1) * (~index);
    end
end

function [z] = upperz_update(LB, UB, x, z, u)
    % Compute the mean x 
    X = sum(x)/size(x,1);
    U = sum(u)/size(u,1);
    res = X+U; 

    % Update the global variable (projection on the box)
    if (res < LB)
        z = LB;
    elseif (res > UB)
        z = UB;
    else
        z = res;
    end
end

function [p] = upper_objective(LB, UB, x, z)
    if (z > UB || z < LB)
        p = NaN;
    else
        p = 0;
    end
end
