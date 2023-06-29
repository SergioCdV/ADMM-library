%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the time-fixed CW optimal L1 problem using ADMM and PVT

close; 
clear; 
clc

%% Define the rendezvous problem and the STM 
% Time span
t = linspace(0,pi,1000); 
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

% Maximum L1 norm
dVmax = 0.2; 

% Cardinality sequence 
K = 50;

%% Define the ADMM problem 
rho = 1;        % AL parameter 

% pre-factor
Atb = pinv(Phi)*b;
pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);

p = repmat(3,1,length(t));
cum_part = cumsum(p);

% Create the functions to be solved 
Obj = @(x,z)(objective(cum_part, z));
X_update = @(x,z,u)(x_update(pInvA, Atb, x, z, u));
Z_update = @(x,z,u)(z_update(dVmax, cum_part, K, rho, x, z, u));

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
dV = reshape(x(:,end), 3, []);
cost = sum(sum(abs(dV),1));
res(:,1) = b-Phi*x(:,end);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval); 
yline(cost, '--')
legend('ADMM cost', 'PVT cost')
grid on;
ylabel('$\Delta V_T$')
xlabel('Iteration $i$')

figure
hold on
yline(dVmax, '--');
stem(t, sum(abs(dV),1), 'filled'); 
grid on;
legend('$\Delta V_{max}$', '$\Delta V_{ADMM}$')
ylabel('$\Delta V$')
xlabel('$t$')

%% Auxiliary functions 
function [x] = x_update(pInvA, Atb, x, z, u) 
   % Impulses update
   x = pInvA*(z-u)+Atb;
end

function [z] = z_update(dVmax, p, K, rho, x, z, u)
    % Impulses update
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:p(i);
        z(sel) = shrinkage(x(sel) + u(sel), 1/rho);     % Minimization
        z(sel) = ball_projection(z(sel), dVmax);        % Ball projection   
        start_ind = p(i) + 1;
    end

    dV = reshape(z,3,[]);
    [~, pos] = sort(sum(abs(dV),1), 'descend');
    
    index = pos(K+1:end);
    for i = 1:length(index)
        z(1 + 3 * (index(i)-1): 3 * index(i)) = zeros(3,1);
    end
end

% Minimization of the L1 norm
function z = shrinkage(x, kappa)
    z = max(0, x-kappa) - max(0, -x-kappa);
end

% Projection onto the L1 ball 
function [x] = ball_projection(x, a)
    if (sum(abs(x)) > a)
        u = sort(x,'descend');
        K = 1:length(x); 

        index = false * ones(1,length(K));
        for i = 1:length(K)
            index(i) = sum(u(1:K(i)))/K(i) < u(K(i));
        end
        index(1) = 1;
        u = u(logical(index));
        rho = sum(u-a)/length(u);
        x = max(x-rho,0);
    end
end

% Objective function
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