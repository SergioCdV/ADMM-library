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
N = 200;

%% Define the rendezvous problem and the STM %%
% Time span
nu = linspace(nu_0, nu_f, N);

% Control input matrix 
B = repmat([0 0; 1 0; 0 0; 0 1], 1, length(nu));

% HCW Phi
Phi = zeros(4, 4 * N);
STM = zeros(4, 4 * N);

Phi1 = [0 2*sin(nu(1)) -3*sin(nu(1)) cos(nu(1)); ...
        0 -2*cos(nu(1)) 3*cos(nu(1)) sin(nu(1)); ...
        0 1 -2 0; ...
        1 3*nu(1) -6*nu(1) 2];

for i = 1:length(nu)
    
    Phi2 = [-2*cos(nu(i)) -2*sin(nu(i)) -3*nu(i) 1; ...
             2*sin(nu(i)) -2*cos(nu(i)) -3 0; ...
             sin(nu(i)) -cos(nu(i)) -2 0; ...
             cos(nu(i)) sin(nu(i)) 0 0];

    Phi(:,1+4*(i-1):4*i) = Phi2;
    STM(:,1+4*(i-1):4*i) = Phi2 * Phi1;
end

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, Phi, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.NeustadtSolver(myMission, myThruster);

rho = 1;                                     % AL parameter 
eps = [1e-5; 1e-4];                                     % Numerical tolerance
[~, sol, ~, myProblem] = myProblem.Solve(eps, rho);

lambda = reshape(sol(1:4), 4, []);
p = reshape(sol(5:end), 2, []);
dV = myProblem.u;

% Optimization
c = (reshape(Phi(:,end-3:end), [4 4])\xf) - (reshape(Phi(:,1:4), [4 4])\x0);
lambda_arz = [0.6481 0 0.6248 -0.0663].';
cost_arz = c.' * lambda_arz;
cost_admm = myProblem.Cost;
Output = myProblem.Report;

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

% Dimensionalization 
% s = s .* repmat([Lc Lc Vc Vc], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * n); 
yline(cost_arz * n, '--')
grid on;
ylabel('$-c^{T} \lambda$')
xlabel('Iteration $i$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
stem(nu, p_norm, 'filled'); 
yline(1, '--')
grid on;
ylabel('$\|\mathbf{p}\|_q$')
xlabel('$\nu$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

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
ylabel('$y$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));