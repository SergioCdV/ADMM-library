%% Optimal Linear Rendezvous via ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Benedikter_2011.m 
% Validated: 

%% HCW, Benedikter 2011 %% 
% Solve for the time-fixed HCW optimal L2/L1 problem using ADMM and PVT

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
mu = 3.986e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [7011e3 0 deg2rad(190) deg2rad(98) 0 0];

nu_0 = 0;            % Initial true anomaly
nu_f = 10;           % Final true anomaly 

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;                 % Initial clock
tf = 7200;              % Final clock

% Initial relative conditions 
x0 = [-pi 1/6 1/4 0];   % In-plane rendezvous
xf = [0 0 0 0];         % Final conditions
   
% Dimensionalization (canonical units)
Lc = Orbit_t(1);        % Characteristic length
Tc = 1/n;               % Characteristic time
Vc = Lc/Tc;             % Characteristic velocity

mu = mu / (Lc^3/Tc^2);  % Gravitational 
n = n * Tc;

x0 = x0.'; 
xf = xf.';

Orbit_t(1) = Orbit_t(1) / Lc;

% Additional parameters 
h = sqrt(mu * Orbit_t(1) * (1-Orbit_t(2)^2));

% Number of possible impulses 
N = 100;

%% Define the rendezvous problem and the STM %%
% Time span
t = linspace(nu_0, nu_f, N);

% Control input matrix 
B = repmat([zeros(2); eye(2)], 1, length(t));

% HCW Phi
Phi = zeros(4, 4 * N);
STM = zeros(4, 4 * N);

Phi1 = [0 -3*sin(t(1)) 2*sin(t(1)) cos(t(1)); ...
        0 -2 1 0; ...
        0  3*cos(t(1)) -2*cos(t(1)) sin(t(1)); ...
        1 -6*t(1) 3*t(1) 2];

for i = 1:length(t)
    
    Phi2 = [-2*cos(t(i)) -3*t(i) -2*sin(t(i)) 1; ...
         sin(t(i)) -2 -cos(t(i)) 0; ...
         2*sin(t(i)) -3 -2*cos(t(i)) 0; ...
         cos(t(i)) 0 sin(t(i)) 0];

    Phi(:,1+4*(i-1):4*i) = Phi2;
    STM(:,1+4*(i-1):4*i) = Phi2 * Phi1;
end

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(t, Phi, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.NeustadtSolver(myMission, myThruster);

iter = 1;
time = zeros(1,iter);

rho = 1/N;                                              % AL parameter 
eps = [1e-5; 1e-4];                                     % Numerical tolerance

for i = 1:iter
    [~, sol, ~, myProblem2] = myProblem.Solve(eps, rho);
    time(i) = myProblem2.SolveTime;
end
myProblem = myProblem2;

lambda = reshape(sol(1:4), 4, []);
p = reshape(sol(5:end), 2, []);
dV = myProblem.u;

%% Outcome

% Norm of the primer vector 
switch (myThruster.p)
    case 'L1'
        dV_norm = sum(abs(dV),1);
    case 'L2'
        dV_norm = sqrt(dot(dV,dV,1));
    case 'Linfty'
        dV_norm = max(abs(dV));
end

switch (myThruster.q)
    case 'L1'
        p_norm = sum(abs(p),1);
    case 'L2'
        p_norm = sqrt(dot(p,p,1));
    case 'Linfty'
        p_norm = max(abs(p));
end

% Impulsive times
ti = dV_norm ~= 0;

% Results
cost_admm = myProblem.Cost;
Output = myProblem.Report;
Time = mean(time);
Nopt = sum(ti);
error = sqrt( dot(myProblem.e, myProblem.e, 1) ); 
t_imp = t(ti);

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(t),4);
s(1,:) = x0.';

% Computation
for i = 1:length(t)
    % Propagate 
    if (i > 1)
        Phi1 = reshape(STM(:,1+4*(i-2):4*(i-1)), [4 4]);
        Phi2 = reshape(STM(:,1+4*(i-1):4*i), [4 4]);
        s(i,:) = s(i-1,:) * (Phi2 * Phi1^(-1)).';
    end

    % Add maneuver
    s(i,3:4) = s(i,3:4) + dV(:,i).';
end

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * Vc); 
grid on;
ylabel('$\mathbf{c}^{T} \mathbf{\lambda}$')
xlabel('Iteration $i$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
scatter(t_imp, ones(1,length(t_imp)), 1e2, 'r', 'Marker', 'x')
legend('$t_i$', 'AutoUpdate', 'off')
plot(t, p_norm, 'b');
yline(1, '--')
grid on;
ylabel('$\|\mathbf{p}\|_q$')
xlabel('$t$')
xlim([t(1) t(end)])

figure
hold on
stem(t, dV_norm * Vc, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_p$ [m/s]')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([t(1) t(end)])

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
hold on
scatter(s(1,1), s(1,2), siz, 'b', 'Marker', 'square');
scatter(s(ti,1), s(ti,2), siz2, 'r', 'Marker', 'x');
scatter(s(end,1), s(end,2), siz, 'b', 'Marker', 'o');
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}_f$', 'AutoUpdate', 'off');
plot(s(:,1), s(:,2), 'b'); 
hold off
xlabel('$x$')
ylabel('$z$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

%% Absolute dynamics figure
S = [s(:,1) -ones(length(t),1)+s(:,2)];

% Rotation
s0 = [-0.83; 0; 0];
o21 = -[0; 0; 1];
for i = 1:length(t)
    o31 = -[cos(t(i)); sin(t(i)); 0];
    o11 = cross(o21, o31);
    Q = [o11 o21 o31];
    aux = Q * [S(i,1); 0; S(i,2)];
    S(i,1:2) = aux([1 2]);
end

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
hold on
scatter(S(1,1), S(1,2), siz, 'b', 'Marker', 'square');
scatter(S(ti,1), S(ti,2), siz2, 'r', 'Marker', 'x');
scatter(S(end,1), S(end,2), siz, 'b', 'Marker', 'o');
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}_f$', 'AutoUpdate', 'off');
plot(S(:,1), S(:,2), 'b'); 
fplot(@(t) cos(t), @(t) sin(t), 'k');
fplot(@(t) cos(t)/1.2, @(t) sin(t)/1.2, 'k--');
hold off
xlabel('$X$')
ylabel('$Y$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));