%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: EUCASS_2023_CW2.m 
% Issue: 0 
% Validated: 

%% Clohessy-Wiltshire rendezvous %% 
% Solve for the time-fixed CW optimal L2 problem using ADMM and PVT

close; 
clear; 
clc

set_graphics();

%% Define the target orbit
% Parameters
mu = 3.986e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [6900e3 0 0 deg2rad(98) 0 0];

% Normalization
a = 1;               % Characteristic length
n = 1;               % Mean motion
Tc = n;              % Characteristic time

%% Define the rendezvous problem and the STM 
tf = 2*pi;                      % Time of flight
    
% Time span
t = linspace(0, tf, 500);

% Relative initial conditions
x0 = [-0.005019; 0.01; 0.01; 0.01; 0.01; 0];  

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

% Residual or missvector
b = xf-M*x0;

%% Optimization
% Define the ADMM problem 
rho = 5e2;        % AL parameter 

% Pre-factor
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
B = -eye(n);        
c = zeros(m,1);

% Problem
Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

% Optimization
tic
[x, z, Output] = Problem.solver();
dV = reshape(x(:,end), 3, []);
cost_admm = sum(sqrt(dot(dV,dV,1)));
toc 

% Pruning
[dV2, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV, 'L2');

%% Outcome 
res(:,1) = b-Phi*x(:,end);
res(:,2) = b-Phi*reshape(dV2, [], 1);
ratio = 1-cost/cost_admm;

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(t),6);
s(1,:) = x0.';

% Computation
for i = 1:length(t)
    % Propagate 
    if (i ~= 1)
        Phi1 = reshape(STM(i-1,:), [6 6])^(-1);
        Phi2 = reshape(STM(i,:), [6 6]);
        s(i,:) = s(i-1,:) * (Phi2 * Phi1).';
    end

    % Add maneuver
    if (norm(dV2(:,i)) ~= 0)
        s(i,4:6) = s(i,4:6) + dV2(:,i).';
    end
end

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval); 
yline(cost, '--')
legend('ADMM cost', 'PVT cost')
grid on;
ylabel('$\Delta V_T$')
xlabel('Iteration $i$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
% yline(dVmax, '--');
stem(t, sqrt(dot(dV,dV,1)), 'filled'); 
grid on;
% legend('$\Delta V_{max}$', '$\Delta V_{ADMM}$')
ylabel('$\Delta V$')
xlabel('$t$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
% yline(dVmax, '--');
stem(t, sqrt(dot(dV2,dV2,1)), 'filled'); 
grid on;
% legend('$\Delta V_{max}$', '$\Delta V_{PVT}$')
ylabel('$\Delta V$')
xlabel('$t$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot3(s(:,1), s(:,2), s(:,3)); 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary functions 
% Proximal minimization for X
function [x] = x_update(pInvA, Atb, x, z, u) 
   x = pInvA*(z-u)+Atb;
end

% Proximal minimization for Z
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

% Proximal minimization of the l2 norm
function z = shrinkage(x, kappa)
    z = max(0, 1 - kappa/norm(x))*x;
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

% Set graphics
function set_graphics()
    %Set graphical properties
   set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end