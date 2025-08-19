%% Optimal Linear Rendezvous via simplex %% 
% Sergio Cuevas del Valle
% Date: 01/06/25
% File: LowThrust.m 
% Issue: 0  

%% Yamanaka-Andersen rendezvous, Arzelier 2016 %% 
% Solve for the time-fixed YA optimal L1 problem using incomplete simplex method %

close; 
clear; 
clc

utils.set_graphics();

%% Define the target orbit and mission parameters
% Parameters
mu = 3.986e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [6763e3 0.0052 0 deg2rad(52) 0 0];

nu_0 = 0;            % Initial true anomaly

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 0;              % Initial clock
tf = 7200;           % Final clocks

K = floor(tf/(2*pi/n));
dt = tf - K * (2*pi/n);

nu_f = 2*pi*K + InverseKeplerEquation(n, Orbit_t(2), nu_0, dt);      % Final true anomaly 

% Initial relative conditions 
x0 = [-30 0.5 8.514e-3 0]*1e3;     % In-plane rendezvous
xf = [-100 0 0 0];                 % Final conditions
   
% Dimensionalization (canonical units)
Lc = Orbit_t(1);        % Characteristic length
Tc = 1/n;               % Characteristic time
Vc = Lc/Tc;             % Characteristic velocity
Ac = Lc/Tc^2;           % Characteristic acceleration

mu = mu / (Lc^3/Tc^2);  % Gravitational 
n = n * Tc;             % Mean motion

x0 = x0 ./ [Lc Lc Vc Vc];
xf = xf ./ [Lc Lc Vc Vc];

x0 = x0.'; 
xf = xf.';

Orbit_t(1) = Orbit_t(1) / Lc;

% Additional parameters 
h = sqrt(mu * Orbit_t(1) * (1-Orbit_t(2)^2));

%% Define the rendezvous problem and the STM %%
% Time span
N = 300;                         % Discretization size
nu = linspace(nu_0, nu_f, N);
t = nu;

K = 0;
for i = 1:length(nu)
    dt = KeplerEquation(n, Orbit_t(2), nu(1), nu(i));
    if (i > 2)
        if (mod(nu(i),2*pi) < mod(nu(i-1),2*pi))
            K = K+1;
        end
    end
    t(i) = 2*K*pi + dt;
end

% Control input matrix 
B = repmat([zeros(2); eye(2)], 1, length(nu));

% YA Phi
L = zeros(4, 4 * N);
Phi = zeros(4, 4 * N);
DT = 0;
K = 0;

for i = 1:length(nu)
    % Constants of motion 
    omega = mu^2 / h^3;                 % True anomaly angular velocity
    k = 1 + Orbit_t(2) * cos(nu(i));    % Transformation
    kp =  - Orbit_t(2) * sin(nu(i));    % Derivative of the transformation

    % Solve Kepler's equation
    dt = KeplerEquation(n, Orbit_t(2), nu(1), nu(i));
    
    % Consider multiple revolutions
    if (i > 2)
        if (mod(nu(i),2*pi) < mod(nu(i-1),2*pi))
            K = K+1;
        end
    else
       Phi0 = YA_Phi(mu, h, Orbit_t(2), DT, nu(1)); 
       invPhi0 = Phi0([1 3 4 6], [1 3 4 6])^(-1);
       L(:,1+4*(i-1):4*i) = [k * eye(2) zeros(2); kp * eye(2) eye(2)/(k * omega)];
    end

    DT = 2*K*pi + dt;
    phi = YA_Phi(mu, h, Orbit_t(2), DT, nu(i));
    
    stm = phi([1 3 4 6], [1 3 4 6]) * invPhi0;

    L(:,1+4*(i-1):4*i) = [k * eye(2) zeros(2); kp * eye(2) eye(2)/(k * omega)];
    STM(:,1+4*(i-1):4*i) = L(:,1+4*(i-1):4*i)^(-1) * stm * L(:,1:4);
end

%% Save parameters for the PS method
save ParametersATV.mat -mat

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, STM, B, x0, xf, K);       % Mission

%% Optimization (until achieving low-thrust) 
dVmin = 0;              % Minimum control authority
dVmax = [Inf 0.2];      % Maximum control authority

dVmax = dVmax / Vc;
dVmin = dVmin / Vc;

Cost = zeros(1,length(dVmax));
Time = zeros(1,length(dVmax));
Error = zeros(1,length(dVmax));
dV = zeros(2 * length(dVmax), N);
dV_norm = zeros(length(dVmax), N);
Max_violation = zeros(1,length(dVmax));

s = zeros(length(nu), 4 * length(dVmax));

for i = 1:length(dVmax)
    % Thruster definition 
    myThruster = thruster( 'L1', dVmin, dVmax(i) );

    % Define the ADMM problem 
    myProblem = RendezvousProblems.GenPotterSolver(myMission, myThruster);

    % Solve the problem
    tic
    [~, dV(1+2*(i-1):2*i,:), ~, myProblem] = myProblem.Solve();
    Time(i) = toc;

    % Outcome 
    switch (myThruster.p)
        case 'L1'
            dV_norm(i,:) = sum(abs(dV(1+2*(i-1):2*i,:)),1);
        case 'L2'
            dV_norm(i,:) = sqrt(dot(dV(1+2*(i-1):2*i,:),dV(1+2*(i-1):2*i,:),1));
        case 'Linfty'
            dV_norm(i,:) = max(abs(dV(1+2*(i-1):2*i,:)));
    end

    switch (myThruster.q)
        case 'L1'
            dV_auth = sum(abs(dV(1+2*(i-1):2*i,:)),1);
        case 'L2'
            dV_auth = sqrt(dot(dV(1+2*(i-1):2*i,:),dV(1+2*(i-1):2*i,:),1));
        case 'Linfty'
            dV_auth = max(abs(dV(1+2*(i-1):2*i,:)));
    end

    % Impulsive times
    ti = dV_norm(i,:) >= 0.01 * max(dV_norm(i,:));
    
    % Results
    Cost(i) = myProblem.Cost;
    Nopt = sum(ti);
    Error(i) = sqrt( dot(myProblem.e, myProblem.e, 1) ); 
    t_imp = t(ti);
    
    cons = dV_auth - dVmax(i);
    cons = cons( abs( cons ) < 0.001 * dVmax(i) );

    if ( ~isempty(cons) )
        Max_violation(i) = max( cons );
    end

    % Preallocation 
    idx = 1+4*(i-1):4*i;

    s(1,idx) = x0.';
    
    % Computation
    for j = 1:length(nu)
        % Propagate 
        if (j > 1)
            Phi1 = reshape(STM(:,1+4*(j-2):4*(j-1)), [4 4]);
            Phi2 = reshape(STM(:,1+4*(j-1):4*j), [4 4]);
            s(j,idx) = s(j-1,idx) * (Phi2 * Phi1^(-1)).';
        end
    
        % Add maneuver
        s(j,idx(end-1:end)) = s(j,idx(end-1:end)) + dV(1+2*(i-1):2*i,j).';
    end
    
    % Dimensionalization 
    s(:,idx) = s(:,idx) .* repmat([Lc Lc Vc Vc], size(s,1), 1);
    s(:,idx) = s(:,idx) / 1e3;
end

save SolutionImpulsiveATV.mat -mat

%% Results 
for i = 1:length(dVmax)
    figure(i)
    hold on
    stem(nu, dV_norm(2,:) * Vc, 'filled', Color=[0 0.4470 0.7410]);
%     stem(nu, dV_norm(2,:) * Vc, 'filled', Color=[0 0.4470 0.7410]); 
end
grid on;
ylabel('$\|\Delta \mathbf{V}\|_1$ [m/s]')
xlabel('$\theta$')
%     legend('$\Delta V_{\mathrm{max}} = \infty\,\mathrm{m/s}$', '$\Delta V_{\mathrm{max}} = 20\,\mathrm{cm/s}$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([nu(1) nu(end)])

siz = repmat(100, 1, 1);
figure 
hold on
scatter(s(1,1), s(1,2), siz, 'b', 'Marker', 'square');
scatter(s(end,1), s(end,2), siz, 'b', 'Marker', 'o');
plot(s(:,1), s(:,2), 'c', 'LineWidth', 1);
plot(s(:,5), s(:,6), 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
legend('$\mathbf{s}_0$', '$\mathbf{s}_f$', '$\Delta V_{\mathrm{max}} = \infty\,\mathrm{m/s}$', '$\Delta V_{\mathrm{max}} = 20\,\mathrm{cm/s}$', 'AutoUpdate', 'off');
hold off
xlabel('$x$ [km]')
ylabel('$z$ [km]')
% xlim([-1.1e3 100])
% ylim([-20 200])
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

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

function [nu_f] = InverseKeplerEquation(n, e, M0, dt)
    % Initial mean anomaly 
    M = M0 + n * dt;

    % Laguerre-Conway's method
    maxIter = 10;       % Maximum number of iterations
    iter = 1;           % Initial iteration
    GoOn = true;        % Convergene boolean flag
    k = 5;              % Laguerre constant
    tol = 1E-15;        % Convergence tolerance

    % Warm start
    u = M + e;
    E = M .* (1-sin(u)) + u .* sin(M) ./ (1+sin(M)-sin(u));

    while (iter < maxIter && GoOn)
        fn = E - e .* sin(E) - M;
        dfn = 1 - e .* cos(E);
        ddfn = e .* sin(E);

        dg(1) = k / (dfn+sqrt( abs((k-1)^2*dfn^2-k*(k-1)*dfn*ddfn)) );
        dg(2) = k / (dfn-sqrt( abs((k-1)^2*dfn^2-k*(k-1)*dfn*ddfn)) );

        dg = dg( abs(dg) == max( abs(dg) ) );

        dn = - fn / dg;
        E = E + dn;

        if all(abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end

    % Final true anomaly 
    sin_nu = (sqrt(1-e^2) * sin(E)) / (1 - e * cos(E));
    cos_nu = (-e + cos(E)) / (1 - e * cos(E));
    nu_f = atan2(sin_nu, cos_nu);
    nu_f = mod(nu_f,2*pi);
end