%% Optimal Linear Rendezvous via simplex %% 
% Sergio Cuevas del Valle
% Date: 03/03/24
% File: SIMBOL.m 
% Issue: 0 
% Validated: 

%% SIMBOL X rendezvous %% 
% Solve for the time-fixed YA optimal using incomplete simplex method %

close; 
clear; 
clc

utils.set_graphics();

%% Define the target orbit and mission parameters
% Parameters
mu = 3.9860e14;       % Gauss constant for the Earth

% Target orbital elements
Orbit_t = [106246.9753e3 0.798788 pi/2 deg2rad(5.2) deg2rad(180) 0];

nu_0 = deg2rad(135);            % Initial true anomaly

% Mean motion
n = sqrt(mu/Orbit_t(1)^3);      

% Mission time
t0 = 7;              % Initial clock
tf = 50002;          % Final clock

% Final true anomaly
K = floor(tf/(2*pi/n));
dt = tf - K * (2*pi/n);

e = Orbit_t(2);
cos_E = (e + cos(nu_0)) / (1 + e * cos(nu_0));
sin_E = (sqrt(1-e^2)*sin(nu_0)) / (1 + e * cos(nu_0));
E = atan2(sin_E, cos_E);
M0 = E - Orbit_t(2)*sin(E);
nu_f = 2*pi*K + InverseKeplerEquation(n, Orbit_t(2), M0, dt);       % Final true anomaly 

% Initial relative conditions 
x0 = [+18309.5 -23764.7 -0.0542 -0.0418];              % In-plane rendezvous
xf = [+335.12 -371.1 +0.00155 +0.0014];                  % Final conditionsxl
   
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
N = 1000;

%% Define the rendezvous problem and the STM %%
% Time span
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

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, STM, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;%0.001 / Vc;                                              % Minimum control authority
dVmax = Inf;%0.3 / Vc;                                            % Maximum control authority
myThruster = thruster('L2', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myProblem = RendezvousProblems.GenPotterSolver(myMission, myThruster);
 
iter = 1; 
time = zeros(1,iter);

for i = 1:iter
    [~, dV, ~, myProblem] = myProblem.Solve();
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
cost = myProblem.Cost;
Time = mean(time);
Nopt = sum(ti);
error = sqrt( dot(myProblem.e, myProblem.e, 1) ); 
nu_imp = nu(ti);
t_imp = t(ti) * Tc;

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
    s(i,3:4) = s(i,3:4) + dV(:,i).';
end

% Reference solution 
[nuref, sref, dV_ref] = ReferenceSolution(Orbit_t, mu, h, n, Vc, x0);

switch (myThruster.p)
    case 'L1'
        dV_norm_ref = sum(abs(dV_ref),1);
    case 'L2'
        dV_norm_ref = sqrt(dot(dV_ref,dV_ref,1));
    case 'Linfty'
        dV_norm_ref = max(abs(dV_ref));
end

% Dimensionalization 
s =    s    .* repmat([Lc Lc Vc Vc], size(s,1), 1);
sref = sref .* repmat([Lc Lc Vc Vc], size(sref,1), 1);

s = s / 1000;
sref = sref / 1000;

%% Results 
figure
hold on
if (dVmax ~= Inf)
    yline(dVmax * Vc * 100, 'k--')
    if (dVmin > 0)
        yline(dVmin, 'k--')
        legend('$\Delta V_{max}$', '$\Delta V_{min}$', 'Autoupdate', 'off')
    else
        legend('$\Delta V_{max}$', 'Autoupdate', 'off')
    end

   switch (myThruster.q)
        case 'Linfty'
            stem(t, max(abs(dV), [], 1) * Vc, 'filled', 'k');
        case 'L1'
            stem(t, sum(abs(dV), [], 1) * Vc, 'filled', 'k');
        otherwise
   end
else
   stem(nuref, dV_norm_ref * Vc * 100, 'filled', 'c'); 
end
stem(nu, dV_norm * Vc * 100, 'filled', Color=[0 0.4470 0.7410]); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_2$ [cm/s]')
xlabel('$\theta$')
legend('Arzelier et al.', 'PS')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([min(nu(1), nuref(1)) max(nu(end), nuref(end))])

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
hold on
scatter(s(1,1), s(1,2), siz, 'b', 'Marker', 'square');
scatter(s(ti,1), s(ti,2), siz2, 'r', 'Marker', 'x');
scatter(s(end,1), s(end,2), siz, 'b', 'Marker', 'o');
plot(sref(:,1), sref(:,2), 'c', 'LineWidth', 0.2); 
plot(s(:,1), s(:,2), 'b', 'LineWidth', 1); 
legend('$\mathbf{s}_0$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}_f$', '$\mathbf{s}_{ref}$', '$\mathbf{s}_{PS}$', 'AutoUpdate', 'off');
hold off
xlabel('$x$ [km]')
ylabel('$z$ [km]')
% xlim([-1.1e3 100])
% ylim([-20 200])
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
stem(nu, abs( (dV_norm - dVmax) ) * Vc * 100 .* (dV_norm > 0 & dV_norm > dVmax), 'filled'); 
stem(nu, abs( (dV_norm - dVmin) ) * Vc * 100 .* (dV_norm > 0 & dV_norm < dVmin), 'filled'); 
grid on;
ylabel('$e$ [cm/s]')
xlabel('$\theta$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([nu(1) nu(end)])

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

    % Newton-method 
    iter = 1; 
    maxIter = 100; 
    GoOn = true; 
    tol = 1e-15;
    E = M;

    while (GoOn && iter < maxIter)
        f = E - e * sin(E) - M; 
        df = 1 - e * cos(E);

        ds = -f/df;
        E = E + ds; 

        if (abs(ds) < tol)
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

function [nu, sref, dV] = ReferenceSolution(Orbit_t, mu, h, n, Vc, x0) 
    % Reference solution Arzelier et al., 2011 and Benedikter, 2019 for L2
    nu_ref = [2.3562 2.7859];                                                    % Impulsive locations
    dV_ref(:,1:2) = [-0.6193 +0.1748; +0.5061 -0.4912] / Vc;              % Non-dimensional impulse sequence

    % Complete the domain
    nu = [];
    dV = [];
    for i = 1:size(nu_ref,2)-1
        aux = linspace(nu_ref(i), nu_ref(i+1), 1000);
        nu = [nu aux];
        dV = [dV dV_ref(:,i) zeros(size(dV_ref,1), size(aux,2)-1)];
        nu = nu(1:end-1);
        dV = dV(:,1:end-1);
    end

    nu = [nu nu_ref(end)];
    dV = [dV dV_ref(:,end)];
    
    % Pre-allocation
    sref = zeros(length(nu),4);
    sref(1,:) = x0.';
    
    t = nu;
    N = length(nu);
    
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
    
    % YA Phi
    L = zeros(4, 4 * N);
    Phi = zeros(4, 4 * N);
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
           Phi0 = YA_Phi(mu, h, Orbit_t(2), 0, nu(1)); 
           invPhi0 = Phi0([1 3 4 6], [1 3 4 6])^(-1);
           L(:,1+4*(i-1):4*i) = [k * eye(2) zeros(2); kp * eye(2) eye(2)/(k * omega)];
        end
    
        DT = 2*K*pi + dt;
        phi = YA_Phi(mu, h, Orbit_t(2), DT, nu(i));
    
        stm = phi([1 3 4 6], [1 3 4 6]) * invPhi0;
        
        L(:,1+4*(i-1):4*i) = [k * eye(2) zeros(2); kp * eye(2) eye(2)/(k * omega)];
        Phi(:,1+4*(i-1):4*i) = L(:,1+4*(i-1):4*i)^(-1) * phi([1 3 4 6], [1 3 4 6]);
        STM(:,1+4*(i-1):4*i) = L(:,1+4*(i-1):4*i)^(-1) * stm * L(:,1:4);
    end

    % Computation
    for i = 1:length(nu)
        % Propagate 
        if (i > 1)
            Phi1 = reshape(STM(:,1+4*(i-2):4*(i-1)), [4 4]);
            Phi2 = reshape(STM(:,1+4*(i-1):4*i), [4 4]);
            sref(i,:) = sref(i-1,:) * (Phi2 * Phi1^(-1)).';
        end
    
        % Add maneuver
        sref(i,3:4) = sref(i,3:4) + dV(:,i).';
    end
end