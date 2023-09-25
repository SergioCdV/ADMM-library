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

Phi = zeros(6, 6 * length(t));
for i = 1:length(t)
    Phi(:,1+6*(i-1):6*i) = reshape(STM(i,:), [6 6]);
end

B = repmat([zeros(3); eye(3)], 1, length(t));

%% Final mission definition 
K = 20;                                                 % Maximum number of impulses
myMission = LinearMission(t, Phi, B, x0, xf, K);        % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = 0.001;                                          % Maximum control authority
myThruster = thruster('L1', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblemPrimal = RendezvousProblems.PrimalSolver(myMission, myThruster);
myProblemHybrid = RendezvousProblems.HybridSolver(myMission, myThruster);

iter = 1;
time = zeros(2,iter);
rho = sqrt(N);                                        % AL parameter 

for i = 1:iter
    [~, dVp, ~, myProblem2p] = myProblemPrimal.Solve(rho);
    time(1,i) = myProblem2p.SolveTime;

    [~, dVh, ~, myProblem2h] = myProblemHybrid.Solve(rho);
    time(2,i) = myProblem2h.SolveTime;
end

dVh = dVh(1:3,:);
myProblemPrimal = myProblem2p;
myProblemHybrid = myProblem2h;

%% Outcome 
switch (myThruster.p)
    case 'L1'
        dV_norm = sum(abs(dVp),1);
        dV2_norm = sum(abs(dVh),1);
    case 'L2'
        dV_norm = sqrt(dot(dVp,dVp,1));
        dV2_norm = sqrt(dot(dVh,dVh,1));
    case 'Linfty'
        dV_norm = max(abs(dVp));
        dV2_norm = max(abs(dVh));
end

% Impulsive times
ti = dV_norm >= 0.01 * max(dV_norm);
ti2 = dV2_norm >= 0.01 * max(dV2_norm);

% Results
cost_h = Vc * sum(dV_norm);
cost_p = Vc * sum(dV2_norm);
Output = myProblemPrimal.Report;
Output2 = myProblemHybrid.Report;

Time =[mean(time(1,:)); mean(time(2,:))];
Nopt = [sum(ti) sum(ti2)];
error(1:2) = sqrt( dot(myProblemPrimal.e, myProblemPrimal.e, 1) ); 
error(3:4) = sqrt( dot(myProblemHybrid.e, myProblemHybrid.e, 1) ); 
t_imp = t(ti);
t_imp2 = t(ti2);

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(t),6);
s(1,:) = x0.';

% Computation
for i = 1:length(t)
    % Propagate 
    if (i > 1)
        Phi1 = reshape(Phi(:,1+6*(i-2):6*(i-1)), [6 6]);
        Phi2 = reshape(Phi(:,1+6*(i-1):6*i), [6 6]);
        s(i,:) = s(i-1,:) * (Phi2 * Phi1^(-1)).';
    end

    % Add maneuver
    s(i,4:6) = s(i,4:6) + dVp(:,i).';
end

% Preallocation 
sh = zeros(length(t),6);
sh(1,:) = x0.';

% Computation
for i = 1:length(t)
    % Propagate 
    if (i > 1)
        Phi1 = reshape(Phi(:,1+6*(i-2):6*(i-1)), [6 6]);
        Phi2 = reshape(Phi(:,1+6*(i-1):6*i), [6 6]);
        sh(i,:) = sh(i-1,:) * (Phi2 * Phi1^(-1)).';
    end

    % Add maneuver
    sh(i,4:6) = sh(i,4:6) + dVh(:,i).';
end

% Dimensionalization 
% s = s .* repmat([Lc Lc Lc Vc Vc Vc], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval * Vc);
plot(1:Output2.Iterations, Output2.objval * Vc);
grid on;
ylabel('$\Delta V_T$ [m/s]')
xlabel('Iteration $i$')
legend('Primal', 'Hybrid')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
yline(dVmax * Vc, '--')
stem(t, dV_norm * Vc * n, 'filled'); 
stem(t, dV2_norm * Vc * n, 'Marker', 'square'); 
grid on;
legend('$\Delta V_{max}$', 'Primal', 'Hybrid')
ylabel('$\|\Delta \mathbf{V}\|_p$ [m/s]')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([0 2*pi])

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
siz3 = repmat(100, sum(ti2), 1);
figure 
view(3)

hold on
scatter3(s(1,1), s(1,2), s(1,3), siz, 'b', 'Marker', 'square');
scatter3(s(ti,1), s(ti,2), s(ti,3), siz2, 'r', 'Marker', 'x');
scatter3(sh(ti2,1), sh(ti2,2), sh(ti2,3), siz3, 'g', 'Marker', '+');
scatter3(s(end,1), s(end,2), s(end,3), siz, 'b', 'Marker', 'o');
plot3(s(:,1), s(:,2), s(:,3), 'b'); 
plot3(sh(:,1), sh(:,2), sh(:,3), 'r'); 
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i^p$', '$\Delta \mathbf{V}_i^h$', '$\mathbf{s}_f$', '$\mathcal{C}_p$', '$\mathcal{C}_h$', 'AutoUpdate', 'off');
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