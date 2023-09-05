%% Optimal Linear Rendezvous via ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Arzelier_2016.m 
% Issue: 0 
% Validated: 

%% Yamanaka-Andersen rendezvous, Arzelier 2016 %% 
% Solve for the time-fixed YA optimal L2/L1 problem using ADMM and PVT

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
nu_f = 8.1832;       % Final true anomaly 

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;              % Initial clock
tf = 7200;           % Final clock

% Initial relative conditions 
x0 = [-30e3 0.5e3 8.514 0];     % In-plane rendezvous
xf = [-100 0 0 0];              % Final conditions
   
% Dimensionalization (canonical units)
Lc = Orbit_t(1);        % Characteristic length
Tc = 1/n;               % Characteristic time
Vc = Lc/Tc;             % Characteristic velocity 

mu = mu / (Lc^3/Tc^2);  % Gravitational 
n = n * Tc;             % Mean motion

x0 = x0 ./ [Lc Lc Vc Vc];
xf = xf ./ [Lc Lc Vc Vc];

x0 = x0.'; 
xf = xf.';

Orbit_t(1) = Orbit_t(1) / Lc;

% Additional parameters 
h = sqrt(mu * Orbit_t(1) * (1-Orbit_t(2)^2));

% Number of possible impulses 
N = 400;

%% Define the rendezvous problem and the STM %%
% Time span
nu = linspace(nu_0, nu_f, N);

% Control input matrix
B = repmat([zeros(2); eye(2)], 1, length(nu));

% YA Phi
L = zeros(4, 4 * N);
Phi = zeros(4, 4 * N);
DT = 0;
k = 0;

for i = 1:length(nu)
    dt = KeplerEquation(n, Orbit_t(2), nu(1), nu(i)); 
    if (i > 2)
        if (mod(nu(i),2*pi) < mod(nu(i-1),2*pi))
            k = k+1;
        end
    end

    DT = 2*k*pi + dt;
    phi = YA_Phi(mu, h, Orbit_t(2), DT, nu(i));

    omega = h / ( (h^2/mu) / (1 + Orbit_t(2)*cos(nu(i))) )^2;
    k = 1 + Orbit_t(2) * cos(nu(i));
    kp = - Orbit_t(2) * sin(nu(i));
    l = [k * eye(3) zeros(3); kp * eye(3) k / omega * eye(3)];
    
    L(:,1+4*(i-1):4*i) = l([1 3 4 6],[1 3 4 6]);
    Phi(:,1+4*(i-1):4*i) = L(:,1+4*(i-1):4*i) \ phi([1 3 4 6], [1 3 4 6]);

    % Control input matrix 
    B(:,1+2*(i-1):2*i) = L(:,1+4*(i-1):4*i) * B(:,1+2*(i-1):2*i);
end

Phi0 = Phi(:,1:4)^(-1);
for i = 1:length(nu)
    STM(:,1+4*(i-1):4*i) = Phi(:,1+4*(i-1):4*i) * Phi0;
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
eps = [1e-3; 1e-4];                                     % Numerical tolerance
[~, sol, ~, myProblem] = myProblem.Solve(eps, rho);

lambda = reshape(sol(1:4), 4, []);
p = reshape(sol(5:end), 2, []);
dV = myProblem.u;

% Optimization
c = (reshape(Phi(:,end-3:end), [4 4])\xf) - (reshape(Phi(:,1:4), [4 4])\x0);
lambda_arz = [-1.177; 1.132; -1.57; 14.36]*1e-4;
cost_arz = -c.' * lambda_arz;
cost_admm = myProblem.Cost * Vc;
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
        s(i,3:4) = s(i,3:4) + dV(:,i).';
    end
end

% Dimensionalization 
s = s .* repmat([Lc Lc Vc Vc], N, 1);

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval); 
yline(cost_arz, '--')
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
stem(nu, dV_norm, 'filled'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|$')
xlabel('$\nu$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot(s(:,1), s(:,2)); 
xlabel('$x$')
ylabel('$y$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary function 
function [dt] = KeplerEquation(n, e, nu_0, nu_f)
    % Initial mean anomaly 
    cos_E = (e + cos(nu_0)) / (1 + e * cos(nu_0));
    sin_E = (sqrt(1-e^2)*sin(nu_0)) / (1 + e * cos(nu_0));
    E = atan2(sin_E, cos_E);
    M0 = E-e*sin(E);
    M0 = mod(M0,2*pi);

    % Final mean anomaly 
    cos_E = (e + cos(nu_f)) / (1 + e * cos(nu_f));
    sin_E = (sqrt(1-e^2)*sin(nu_f)) / (1 + e * cos(nu_f));
    E = atan2(sin_E, cos_E);
    Mf = E-e*sin(E);
    Mf = mod(Mf,2*pi);

    % Time step 
    dM = Mf-M0;
    if (dM < 0)
        Mf = Mf + 2 * pi;
        dM = Mf-M0;
    end
    dt = dM/n;
end