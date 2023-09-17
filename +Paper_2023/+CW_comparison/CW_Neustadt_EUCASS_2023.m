%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 10/09/23
% File: EUCASS_2023_CW.m 
% Issue: 0 
% Validated: 

%% Clohessy-Wiltshire rendezvous %% 
% Solve for the time-fixed CW optimal Lp problem using ADMM and PVT

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
Lc = Orbit_t(1);     % Charactertisc distance
n = sqrt(mu/Lc^3);   % Mean motion
Tc = 1/n;            % Characteristic time
Vc = Lc/Tc;          % Characteristic velocity
n = 1;               % Characteristic frequency

%% Define the rendezvous problem and the STM 
tf = 2*pi;                      % Time of flight
    
% Time span
N = 250;
t = linspace(0, tf, N);
nu = t;

% Relative initial conditions
x0 = [-0.005019; 0.01; 0.01; 0.01; 0.01; 0]; 

% Relative final conditions 
xf = zeros(6,1);

% CW STM
STM = zeros(6, 6 * N);

CW_STM(:,1:6) = [4-3*cos(n*t.') -6*n*t.'+6*sin(n*t.') zeros(length(t),1) 3*n*sin(n*t.') -6*n+6*n*cos(n*t.') zeros(length(t),1)]; 
CW_STM(:,7:12) = [zeros(length(t),1) ones(length(t),1) zeros(length(t),4)];
CW_STM(:,13:18) = [zeros(length(t),2) cos(n*t.') zeros(length(t),2) -n*sin(n*t.')];
CW_STM(:,19:24) = [sin(n*t.')/n -2/n+2*cos(n*t.')/n zeros(length(t),1) cos(n*t.') -2*sin(n*t.') zeros(length(t),1)];
CW_STM(:,25:30) = [2/n-2*cos(n*t.')/n 4/n*sin(n*t.')-3*t.' zeros(length(t),1) 2*sin(n*t.') -3+4*cos(n*t.') zeros(length(t),1)];
CW_STM(:,31:36) = [zeros(length(t),2) sin(n*t.')/n zeros(length(t),2) cos(n*t.')];

Phi = zeros(6, 6 * length(t));

for i = 1:length(nu)
    
    Phi2 = [-2*cos(nu(i)) -2*sin(nu(i)) -3*nu(i) 1; ...
             2*sin(nu(i)) -2*cos(nu(i)) -3 0; ...
             sin(nu(i)) -cos(nu(i)) -2 0; ...
             cos(nu(i))  sin(nu(i)) 0 0];

    Phi3 = [cos(nu(i)) sin(nu(i)); -sin(nu(i)) cos(nu(i))];

    Phi(:,1+6*(i-1):6*i) = reshape(CW_STM(i,:), [6 6]);%blkdiag(Phi2, Phi3 * Phi4 );
    STM(:,1+6*(i-1):6*i) = reshape(CW_STM(i,:), [6 6]);%blkdiag(Phi2 * Phi1, Phi3 * (Phi4\eye(2)) );
end

B = repmat([zeros(3); eye(3)], 1, length(t));

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, Phi, B, x0, xf, K);        % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('Linfty', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.NeustadtSolver(myMission, myThruster);

iter = 1;
time = zeros(1,iter);
rho = 1/N;                                              % AL parameter 
eps = [1e-6; 1e-5];                                     % Numerical tolerance

for i = 1:iter
    [~, sol, ~, myProblem2] = myProblem.Solve(eps, rho);
    time(i) = myProblem2.SolveTime;
end
myProblem = myProblem2;

lambda = reshape(sol(1:6), 6, []);
p = reshape(sol(7:end), 3, []);
dV = myProblem.u;

% Pruning
[dV2, cost] = PVT_pruner(STM, [zeros(3); eye(3)], dV, 'L2');

%% Outcome 
switch (myThruster.p)
    case 'L1'
        dV_norm = sum(abs(dV),1);
        dV2_norm = sum(abs(dV2),1);
    case 'L2'
        dV_norm = sqrt(dot(dV,dV,1));
        dV2_norm = sqrt(dot(dV2,dV2,1));
    case 'Linfty'
        dV_norm = max(abs(dV));
        dV2_norm = max(abs(dV2));
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
cost_pvt = sum(dV2_norm) * Vc;
cost_admm = myProblem.Cost * Vc;
Output = myProblem.Report;
Time = mean(time);
Nopt = sum(ti);
error = sqrt( dot(myProblem.e, myProblem.e, 1) ); 
t_imp = t(ti);

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(t),6);
s(1,:) = x0.';

% Computation
for i = 1:length(t)
    % Propagate 
    if (i > 1)
        Phi1 = reshape(STM(:,1+6*(i-2):6*(i-1)), [6 6]);
        Phi2 = reshape(STM(:,1+6*(i-1):6*i), [6 6]);
        s(i,:) = s(i-1,:) * (Phi2 * Phi1^(-1)).';
    end

    % Add maneuver
    s(i,4:6) = s(i,4:6) + dV(:,i).';
end

% Dimensionalization 
% s = s .* repmat([Lc Lc Lc Vc * n Vc * n Vc * n], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * Vc * n); 
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
stem(t, dV_norm * Vc * n, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_p$ [m/s]')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([0 2*pi])

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
view(3)
hold on
scatter3(s(1,1), s(1,2), s(1,3), siz, 'b', 'Marker', 'square');
scatter3(s(ti,1), s(ti,2), s(ti,3), siz2, 'r', 'Marker', 'x');
scatter3(s(end,1), s(end,2), s(end,3), siz, 'b', 'Marker', 'o');
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}_f$', 'AutoUpdate', 'off');
plot3(s(:,1), s(:,2), s(:,3), 'b'); 
hold off
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary functions
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