%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the time-fixed CW optimal L2 problem using ADMM and PVT

close; 
clear; 
clc

%% Define the rendezvous problem and the STM 
TOF = 2 * pi;   % Time of flight

%% Optimization
% Time span
t = linspace(0,TOF,5);
n = 1;

% Relative initial conditions 
x0 = [1; 2; 3; -1; -pi; 2]; 

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
    Phi(:,1+3*(i-1):3*i) = reshape(STM(i,:), [6 6])^(-1) * [zeros(3); eye(3)];
end

% Residual or missvector
b = [M^(-1)*xf-x0; zeros(3 * length(t),1)];

%% Define the ADMM problem 
rho = 1;        % AL parameter 
    
% Constraint 
[m, n] = size(Phi); 
A = eye(m + 3 * length(t));
B = -eye(m + 3 * length(t));        
c = zeros(m + 3 * length(t),1);

% Create the functions to be solved 
Obj = @(x,z)(objective(b, x));
X_update = @(x,z,u)(x_update(b, rho, x, z, u));
Z_update = @(x,z,u)(z_update(Phi, rho, x, z, u));

% Problem
Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
Problem.alpha = 1;

%% Optimization
[x, z, Output] = Problem.solver();
lambda = x(1:6,end);
p = reshape(x(7:end,end), 3, []);

% for i = 1:length(t)
%     p(:,i) = Phi(:,1+3*(i-1):3*i).' * lambda;
% end

np = sqrt(dot(p,p,1));
plot(t,np)

%% 
% Get the impulses and the locations
G = zeros(6,6);
for i = 1:length(t)
    D = M * reshape(STM(i,:), [6 6])^(-1) * [zeros(3); eye(3)];
    G = G + norm(dV(:,i)) * (D*D.');
end

lambda = G\b;

mag_p = sqrt(dot(p,p,1));
error = abs(1-mag_p);

dV = reshape(x(1:n,end), 3, []);
cost_admm = sum(sqrt(dot(dV,dV,1)));

%% Pruning
[dV2, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV, 'L2');

%% Outcome 
res(:,1) = b-Phi*x(1:n,end);
res(:,2) = b-Preshape(dV2, [], 1);
ratio = 1-cost/cost_admm;

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
stem(t, sqrt(dot(dV,dV,1)), 'filled'); 
grid on;
legend('$\Delta V_{max}$', '$\Delta V_{ADMM}$')
ylabel('$\Delta V$')
xlabel('$t$')

figure
hold on
yline(dVmax, '--');
stem(t, sqrt(dot(dV2,dV2,1)), 'filled'); 
grid on;
legend('$\Delta V_{max}$', '$\Delta V_{PVT}$')
ylabel('$\Delta V$')
xlabel('$t$')

figure
hold on
plot(t, mag_p); 
grid on;
ylabel('$\|p\|$')
xlabel('$t$')

figure
hold on
plot(t, log(error)); 
grid on;
ylabel('$e_p$')
xlabel('$t$')

if (parametric)
    figure
    plot(TOF, cost);
    grid on;
    ylabel('cost')
    xlabel('TOF')
end

%% Auxiliary functions 
function [x] = x_update(c, rho, x, z, u) 
   % Lagrange multiplier update
   x = (z-u)+c/rho;
end

function [z] = z_update(Phi, rho, x, z, u)
    % Lagrange update (projection on the PSD cone)
    A = zeros(size(Phi,2),6);
    b = zeros(size(Phi,2),1);
    p = zeros(size(Phi,2),1);

    for i = 1:size(Phi,2)/3
        % Update of the primer vector 
        p(1+3*(i-1):3*i) = Phi(:,1+3*(i-1):3*i).' * (x(1:6,1)+u(1:6,1));
        pn = norm( p(1+3*(i-1):3*i) );
        if (pn > 1)
            p(1+3*(i-1):3*i) = p(1+3*(i-1):3*i)/pn;
        end
            
        % Linear system
        b(1+3*(i-1):3*i,:) = p(1+3*(i-1):3*i);
        A(1+3*(i-1):3*i,:) = Phi(:,1+3*(i-1):3*i).';
    end

    % Z-update
    pInvA = (eye(size(A,2))-pinv(A)*A); 
    z(1:6,1) = pInvA*( x(1:6)+u(1:6) )+A\b;
    z(7:end) = p;
end

% Euclidean norm proximal operator
function z = shrinkage(x, kappa)
    % L2 minimization
    z = max(0, 1 - kappa/norm(x))*x;
end

% Objective
function [p] = objective(c, x)
    p = -dot(c,x);
end