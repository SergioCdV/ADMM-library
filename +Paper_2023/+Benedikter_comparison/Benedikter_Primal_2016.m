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
Orbit_t = [6763e3 0.0052 0 deg2rad(52) 0 0];

nu_0 = 0;            % Initial true anomaly
nu_f = 10;           % Final true anomaly 

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;                 % Initial clock
tf = 7200;              % Final clock

% Initial relative conditions 
x0 = [-pi 1/4 1/6 0];   % In-plane rendezvous
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
nu = linspace(nu_0, nu_f, N);

% Control input matrix 
B = repmat([0 0; 1 0; 0 0; 0 1], 1, length(nu));

% HCW Phi
STM = zeros(4, 4 * N);

i = 1;
Phi1 = [0 2*sin(nu(i)) -3*sin(nu(i)) cos(nu(i)); ...
        0 -2*cos(nu(i)) 3*cos(nu(i)) sin(nu(i)); ...
        0 1 -2 0; ...
        1 3*nu(i) -6*nu(i) 2];

for i = 1:length(nu)
    
    Phi2 = [-2*cos(nu(i)) -2*sin(nu(i)) -3*nu(i) 1; ...
             2*sin(nu(i)) -2*cos(nu(i)) -3 0; ...
             sin(nu(i)) -cos(nu(i)) -2 0; ...
             cos(nu(i))  sin(nu(i)) 0 0];

    STM(:,1+4*(i-1):4*i) = Phi2 * Phi1;

    pSTM(i,:) = reshape(STM(:,1+4*(i-1):4*i), 1, []);
end

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, STM, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.PrimalSolver(myMission, myThruster);

rho = 1;    % AL parameter 

[~, dV, ~, myProblem] = myProblem.Solve(rho);

% Optimization
cost_admm = myProblem.Cost * n;
Output = myProblem.Report;

% Pruning 
[dV2, cost] = PVT_pruner(pSTM, B(:,1:2), dV, 'L2');

%% Outcome 
switch (myThruster.p)
    case 'L1'
        dV_norm = sum(abs(dV),1);
    case 'L2'
        dV_norm = sqrt(dot(dV,dV,1));
    case 'Linfty'
        dV_norm = max(abs(dV));
end

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(nu),4);
s(1,:) = x0.';

% Computation
for i = 1:length(nu)
    % Propagate 
    if (i > 1)
        Phi1 = reshape(STM(:,1+4*(i-2):4*(i-1)), [4 4]);
        Phi2 = reshape(STM(:,1+4*(i-1):4*i), [4 4]);
        s(i,:) = s(i-1,:) * (Phi2 * Phi1^(-1)).';
    end

    % Add maneuver
    if (norm(dV(:,i)) ~= 0)
        s(i,[2 4]) = s(i,[2 4]) + dV(:,i).';
    end
end

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * n); 
grid on;
ylabel('$\Delta V_T$')
xlabel('Iteration $i$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
stem(nu, dV_norm * n, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|$')
xlabel('$\nu$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot(s(:,1), s(:,3)); 
xlabel('$x$')
ylabel('$z$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));