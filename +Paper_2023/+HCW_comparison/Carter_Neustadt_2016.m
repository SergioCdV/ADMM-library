%% Optimal Linear Rendezvous via ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Kara_2011.m 
% Validated: 

%% HCW, Example I, Kara Zaitri 2011 %% 
% Solve for the time-fixed HCW optimal L2/L1 problem using ADMM and PVT

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
mu = 3.986e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [6763e3 0 deg2rad(0) deg2rad(0) 0 0];

nu_0 = 0;            % Initial true anomaly
nu_f = 2*pi;         % Final true anomaly 

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;                 % Initial clock
tf = 7200;              % Final clock

% Initial relative conditions 
x0 = [1 0 0 0];         % In-plane rendezvous
xf = [0 0 0 0.427];     % Final conditions
   
% Dimensionalization (canonical units)
Lc = Orbit_t(1);        % Characteristic length
Tc = 1/n;               % Characteristic time
Vc = Lc/Tc;             % Characteristic velocity

mu = mu / (Lc^3/Tc^2);  % Gravitational 

% x0 = x0 ./ [Lc Lc Vc Vc];
% xf = xf ./ [Lc Lc Vc Vc];

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
B = repmat([0 0; 1 0; 0 0; 0 1], 1, length(t));

% HCW Phi
Phi = zeros(4, 4 * N);
STM = zeros(4, 4 * N);

Phi1 = [0 2*sin(t(1)) -3*sin(t(1)) cos(t(1)); ...
        0 -2*cos(t(1)) 3*cos(t(1)) sin(t(1)); ...
        0 1 -2 0; ...
        1 3*t(1) -6*t(1) 2];

for i = 1:length(t)
    
    Phi2 = [-2*cos(t(i)) -2*sin(t(i)) -3*t(i) 1; ...
             2*sin(t(i)) -2*cos(t(i)) -3 0; ...
             sin(t(i)) -cos(t(i)) -2 0; ...
             cos(t(i)) sin(t(i)) 0 0];

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
c = (reshape(Phi(:,end-3:end), [4 4])\xf) - (reshape(Phi(:,1:4), [4 4])\x0);
lambda_arz = [0.6481 0 0.6248 -0.0663].';
cost_arz = c.' * lambda_arz;

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
    s(i,[2 4]) = s(i,[2 4]) + dV(:,i).';
end

% Dimensionalization 
% s = s .* repmat([Lc Lc Vc Vc], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval); 
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
xlim([0 2*pi])

figure
hold on
stem(t, dV_norm, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_p$ [m/s]')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([0 2*pi])

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
hold on
scatter(s(1,1), s(1,3), siz, 'b', 'Marker', 'square');
scatter(s(ti,1), s(ti,3), siz2, 'r', 'Marker', 'x');
scatter(s(end,1), s(end,3), siz, 'b', 'Marker', 'o');
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}_f$', 'AutoUpdate', 'off');
plot(s(:,1), s(:,3), 'b'); 
hold off
xlabel('$x$')
ylabel('$y$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));