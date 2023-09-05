%% Optimal Linear Rendezvous via ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Arzelier_2016.m 
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
tf = 1.5*pi;                      % Time of flight
    
% Time span
t = linspace(0, tf, 5e2);

% Relative initial conditions
x0 = [-0.005019; 0.01; 0.01; 0.01; 0.01; 0];  
% x0 = [1; 2; 3; -1; -pi; 2];

% Relative final conditions 
xf = zeros(6,1);

% CW STM
STM(:,1:6) = [4-3*cos(n*t.') -6*n*t.'+6*sin(n*t.') zeros(length(t),1) 3*n*sin(n*t.') -6*n+6*n*cos(n*t.') zeros(length(t),1)]; 
STM(:,7:12) = [zeros(length(t),1) ones(length(t),1) zeros(length(t),4)];
STM(:,13:18) = [zeros(length(t),2) cos(n*t.') zeros(length(t),2) -n*sin(n*t.')];
STM(:,19:24) = [sin(n*t.')/n -2/n+2*cos(n*t.')/n zeros(length(t),1) cos(n*t.') -2*sin(n*t.') zeros(length(t),1)];
STM(:,25:30) = [2/n-2*cos(n*t.')/n 4/n*sin(n*t.')-3*t.' zeros(length(t),1) 2*sin(n*t.') -3+4*cos(n*t.') zeros(length(t),1)];
STM(:,31:36) = [zeros(length(t),2) sin(n*t.')/n zeros(length(t),2) cos(n*t.')];

Phi = zeros(6, 6 * length(t));
for i = 1:length(t)
    Phi(:,1+6*(i-1):6*i) = reshape(STM(i,:), [6 6]);
end

% Control input matrix 
B = repmat([zeros(3); eye(3)], 1, length(t));

%% Final mission definition 
K = 6;                                                % Maximum number of impulses
myMission = LinearMission(t, Phi, B, x0, xf, K);        % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.NeustadtSolver(myMission, myThruster);

rho = 1/length(t);                                      % AL parameter 
eps = [1e-3; 1e-4];                                     % Numerical tolerance
[~, sol, ~, myProblem] = myProblem.Solve(eps, rho);

p = reshape(sol(7:end), 3, []);
dV = myProblem.u;

% Optimization
cost_admm = myProblem.Cost;
Output = myProblem.Report;

% Pruning
[dV2, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV, 'L1');

%% Outcome 
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
        s(i,4:6) = s(i,4:6) + dV(:,i).';
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
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
% yline(dVmax, '--');
stem(t, sqrt(dot(dV2,dV2,1)), 'filled'); 
grid on;
% legend('$\Delta V_{max}$', '$\Delta V_{PVT}$')
ylabel('$\Delta V$')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot3(s(:,1), s(:,2), s(:,3)); 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));