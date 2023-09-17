%% Optimal Linear Rendezvous via ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Kara_2011.m 
% Validated: 

%% HCW, Kara Zaitri 2011 %% 
% Solve for the time-fixed HCW optimal L2/L1 problem using ADMM and PVT

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
mu = 3.986e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [6763e3 0.0052 0 deg2rad(52) 0 0];

nu_0 = 0;            % Initial true anomaly
nu_f = 24*pi;        % Final true anomaly 

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;                 % Initial clock
tf = nu_f/n;            % Final clock

% Initial relative conditions 
x0 = [-10e3 0 0 0];     % In-plane rendezvous
xf = [-100 0 0 0];      % Final conditions
   
% Dimensionalization (canonical units)
Lc = Orbit_t(1);        % Characteristic length
Tc = 1/n;               % Characteristic time
Vc = Lc/Tc;             % Characteristic velocity 

mu = mu / (Lc^3/Tc^2);  % Gravitational 

x0 = x0 ./ [Lc Lc Vc Vc];
xf = xf ./ [Lc Lc Vc Vc];

x0 = x0.'; 
xf = xf.';

Orbit_t(1) = Orbit_t(1) / Lc;

% Additional parameters 
h = sqrt(mu * Orbit_t(1) * (1-Orbit_t(2)^2));

% Number of possible impulses 
N = 300;

%% Define the rendezvous problem and the STM %%
% Time span
t = linspace(nu_0, nu_f, N);

% Control input matrix 
B = repmat([0 0; 1 0; 0 0; 0 1], 1, length(t));

% HCW Phi
STM = zeros(4, 4 * N);

i = 1;
Phi1 = [0 2*sin(t(i)) -3*sin(t(i)) cos(t(i)); ...
        0 -2*cos(t(i)) 3*cos(t(i)) sin(t(i)); ...
        0 1 -2 0; ...
        1 3*t(i) -6*t(i) 2];

for i = 1:length(t)
    
    Phi2 = [-2*cos(t(i)) -2*sin(t(i)) -3*t(i) 1; ...
             2*sin(t(i)) -2*cos(t(i)) -3 0; ...
             sin(t(i)) -cos(t(i)) -2 0; ...
             cos(t(i))  sin(t(i)) 0 0];

    STM(:,1+4*(i-1):4*i) = Phi2 * Phi1;
end

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(t, STM, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.PrimalSolver(myMission, myThruster);

rho = N;    % AL parameter 

iter = 1; 
time = zeros(1,iter);

for i = 1:iter
    [~, dV, ~, myProblem] = myProblem.Solve(rho);
    time(i) = myProblem.SolveTime;
end

%% Outcome 
switch (myThruster.p)
    case 'L1'
        dV_norm = sum(abs(dV),1);
    case 'L2'
        dV_norm = sqrt(dot(dV,dV,1));
    case 'Linfty'
        dV_norm = max(abs(dV));
end

% Impulsive times
ti = dV_norm >= 0.01 * max(dV_norm);

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
    s(i,[2 4]) = s(i,[2 4]) + dV(:,i).';
end

% Dimensionalization 
% s = s .* repmat([Lc Lc Vc Vc], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * Vc * 100); 
grid on;
ylabel('$\Delta V_T$ [cm/s]')
xlabel('Iteration $i$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
stem(t, dV_norm * Vc * 100, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_p$ [cm/s]')
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