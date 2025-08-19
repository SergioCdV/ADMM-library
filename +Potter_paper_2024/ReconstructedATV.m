%% Optimal Linear Rendezvous via simplex %% 
% Sergio Cuevas del Valle
% Date: 03/06/25
% File: LowThrust.m 
% Issue: 0  

%% Yamanaka-Andersen rendezvous, Arzelier 2016 %% 
% Solve for the time-fixed YA optimal L1 problem using incomplete simplex method %

close; 
clear; 
clc

utils.set_graphics();

%% Initial data 
load('RendezvousATVSolution.mat')

% Low-thrust solution
Xlt = x; 
Ult = u; 
nu = tau(1,:);

% Regression of the low-thrust trajectory 
order = 20; 
pol = zeros(size(Xlt,1), order+1);
for i = 1:size(Xlt,1)
    pol(i,:) = polyfit(nu, Xlt(i,:), order);
end

%% Define the rendezvous problem and the STM %%
N = 300; 
nu = linspace(nu_0, nu_f, N);
t = nu;

STM = zeros(4, 4 * N);
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

%% Optimization (equivalent impulsive)
dVmin = 0;        % Minimum control authority
dVmax = [3 7];  % Maximum control authority

dVmax = dVmax / Vc;
dVmin = dVmin / Vc;

% Initial solution 
m = size(Xlt,1); 
n = size(Ult,1);
b = zeros(m * (N - 1), 1); 
A = zeros(m * (N - 1), N * n); 

sampled_x = zeros(m, N); 

for i = 1:m
    sampled_x(i,:) = polyval(pol(i,:), nu);
end

sampled_x = sampled_x * 1000;
for j = 1:N
    sampled_x(:,j) = sampled_x(:,j) ./ [Lc; Lc; Vc; Vc];
end

for i = 2:N
%     if ( i == N )
%         A(1 + m * (i-2) : m * (i-1), 1 + n * (i-1) : n * i) = B(:,1:n);
%     end

    Phi2 = reshape(STM(:,1+m*(i-1):m*i), [m m]); 
    
    for j = 1:i-1
        Phi1 = reshape(STM(:,1+m*(j-1):m*j), [m m]);
        Phi = Phi2 * Phi1^(-1);
        A(1 + m * (i-2) : m * (i-1), 1 + n * (j-1) : n * j) = Phi * B(:,1:n);

        if ( j == 1 )
            b(1 + m * (i-2) : m * (i-1), 1) = sampled_x(:,i) - Phi * sampled_x(:,1);
        end
    end
end

dV0 = A \ b; 
dV0 = reshape(dV0, n, []);

s = zeros(N, m * length(dVmax));
Time = zeros(1, length(dVmax)); 
Cost = zeros(1, length(dVmax));
dV_norm = zeros(length(dVmax), N);

for i = 1:length(dVmax)
    % Thruster definition 
    myThruster = thruster( 'L1', dVmin, dVmax(i) );

    % Define the ADMM problem 
    myProblem = RendezvousProblems.GenPotterSolver(myMission, myThruster);
        
    % Solve the problem
    tic
    [dV, cost] = myProblem.PVT_pruner(STM, B, dV0, myThruster.umax, myThruster.umin, myThruster.p, false);
    Time(i) = toc;

    % Outcome 
    switch (myThruster.p)
        case 'L1'
            dV_norm(i,:) = sum(abs(dV),1);
        case 'L2'
            dV_norm(i,:) = sqrt(dot(dV,dV,1));
        case 'Linfty'
            dV_norm(i,:) = max(abs(dV));
    end

    switch (myThruster.q)
        case 'L1'
            dV_auth = sum(abs(dV),1);
        case 'L2'
            dV_auth = sqrt(dot(dV,dV,1));
        case 'Linfty'
            dV_auth = max(abs(dV));
    end

    cons = dV_auth - dVmax(i);
    cons = cons( abs( cons ) < 0.001 * dVmax(i) );
    
    if ( ~isempty(cons) )
        Max_violation = max( cons );
    end
    
    % Results
    Cost(i) = sum(dV_norm(i,:)) * Vc;

    % Trajectory
    idx = 1 + m * (i-1) : m * i;

    s(1,idx) = Xlt(:,1).' ./ [Lc Lc Vc Vc] * 1000;
    
    for j = 1:length(nu)
        % Propagate 
        if (j > 1)
            Phi1 = reshape(STM(:,1+4*(j-2):4*(j-1)), [4 4]);
            Phi2 = reshape(STM(:,1+4*(j-1):4*j), [4 4]);
            s(j,idx) = s(j-1,idx) * (Phi2 * Phi1^(-1)).';
        end
    
        % Add maneuver
        s(j,idx(end-1:end)) = s(j,idx(end-1:end)) + dV(:,j).';
    end

    % Dimensionalization 
    s(:,idx) = s(:,idx) .* repmat([Lc Lc Vc Vc], size(s,1), 1);
    s(:,idx) = s(:,idx) / 1e3;
end

sampled_x = sampled_x / 1000;
for j = 1:N
    sampled_x(:,j) = sampled_x(:,j) .* [Lc; Lc; Vc; Vc];
end

%% Results 
siz = repmat(100, 1, 1);
figure 
hold on
scatter(s(1,1), s(1,2), siz, 'b', 'Marker', 'square');
scatter(s(end,1), s(end,2), siz, 'b', 'Marker', 'o');
plot(sampled_x(1,:), sampled_x(2,:), 'c', 'LineWidth', 1);
plot(s(:,1), s(:,2), 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
plot(s(:,5), s(:,6), 'k', 'LineWidth', 1);
legend('$\mathbf{s}_0$', '$\mathbf{s}_f$', 'LT', '$\Delta V_{\mathrm{max}} = 3\,\mathrm{m/s}$', '$\Delta V_{\mathrm{max}} = 7\,\mathrm{m/s}$', 'AutoUpdate', 'off');
hold off
xlabel('$x$ [km]')
ylabel('$z$ [km]')
% xlim([-1.1e3 100])
% ylim([-20 200])
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

%%
figure
hold on
stem(nu, dV_norm(1,:) * Vc, 'filled', Color=[0 0.4470 0.7410]);
stem(nu, dV_norm(2,:) * Vc, 'filled', 'k');
legend('$\Delta V_{\mathrm{max}} = 3\,\mathrm{m/s}$', '$\Delta V_{\mathrm{max}} = 7\,\mathrm{m/s}$')
grid on;
ylabel('$\|\Delta \mathbf{V}\|_1$ [m/s]')
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