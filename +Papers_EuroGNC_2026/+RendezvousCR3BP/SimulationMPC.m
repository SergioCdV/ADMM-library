%% Optimal Linear Rendezvous via ADMM in the CR3BP %% 
% Sergio Cuevas del Valle
% Date: 14/10/25
% File: SimulationMPC.m 
% Issue: 0 
% Validated: 

%% C3BP rendezvous via MPC % 
% Solve for a time-fixed optimal rendezvous in the CR3BP (EM L2) via ADMM, MPC and PVT %

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
Lc = 384399e3;          % Characteristic length
Tc = 2.361e6;           % Characteristic time
mu = 0.0121505856;      % Gravitational constan of the system
c2 = 3.190425213622208; % Richardson coefficient
c2 = 4.284390615378499;

% Dimensionalization (canonical units)
Vc = Lc / Tc * 2*pi;    % Characteristic velocity 

% Mission time
nu_0 = 0;           % Initial clock
nu_f = 2.75;           % Final clock

% Initial relative conditions 
x0 = [0.824024728136525; 0; -0.054501847320725; 0; 0.164671964079122; 0];          

% Final conditions
xf = [0.823639438925721; 0; +0.043281569720089; 0; 0.152567980892620; 0];

% Co-orbital initial conditions
x0 = [x0; xf - x0];
xf = [xf; zeros(6,1)];

% Number of possible impulses 
N = 100;                        % Number of steps
Ts = (nu_f - nu_0) / N;

%% Mission definition
% Maximum number of impulses
K = Inf;                                                       

% Control input matrix
B = [zeros(3); eye(3)];

%% Thruster definition 
dVmin = 0;       % Minimum control authority
dVmax = Inf;     % Maximum control authority

myThruster = Actuator(src.VectorNorm.L2, dVmin, dVmax);

%% Optimization
% myPrimalProblem = MinimumNormSolvers.PrimalSolver(myMission, myThruster);

% Optimization
eps = [1e-6; 1e-5];                                         % Numerical tolerance
options = odeset('RelTol', 2.25E-14, 'AbsTol', 1E-22);      % Numerical integration

dV_final = [];
S = [];

while ( N > 0 )
    % Update of the mission
    nu = linspace(nu_0, nu_f, N);
    Phi = wrapper_STM(c2, nu_0, nu_f, N);
    Binp = repmat( B, 1, N );
    myMission = Missions.FuelMission(nu, Phi, Binp, x0(7:end,1), xf(7:end,1), K);
    
    % Pre-allocation of the impulses
    dV = zeros(3,N);                   

    % Define the ADMM problem 
    rho = 1 / N;                       % AL parameter 
    myDualProblem = MinimumNormSolvers.NeustadtSolver(myMission, myThruster);
    
    % Optimization
    [~, sol, ~, myDualProblemSolved] = myDualProblem.Solve(eps, rho^(3/2));
    dV(1:3,:) = myDualProblemSolved.u;

    % Primal resolution
%     [~, dV(4:6,:), ~, myPrimalProblemSolved] = myPrimalProblem.Solve( 1/rho );
    
    % Apply the first impulse 
    x0(10:12,1) = x0(10:12,1) + dV(:,1);
    S = [S x0];
    dV_final = [dV_final dV(:,1)];

    % Integration the coasting solution
    tspan = [nu_0 nu_0 + Ts];
    [~, s] = ode113( @(t,s)cr3bp_equations(mu, t, s), tspan, x0, options );
    
    % Update the IVP
    nu_0 = nu_0 + Ts;           % Update the initial clock
    x0 = s(end,:).';            % New relative initial conditions

    % Update of the MPC loop 
%     N = N - 1;                  % Reduce the number of optimization steps

    disp(N)
end

%% Outcome
% Cost function
dV_norm = myThruster.p.ComputeVectorNorm( dV_final(1:3,:) );

% Impulsive times 
ti = dV_norm ~= 0;

% Results
cost = sum( dV_norm(1,ti), 2 ) * Vc;
% error = sqrt( dot(myProblemSolved.e, myProblemSolved.e, 1) );
Nopt = sum(ti,2);

%% Save results 
% save +Papers_EuroGNC_2026\+RendezvousCR3BP\MPC_L2

%% Dimensionalizations 
N = size(S,2);
dim = [Lc Lc Lc Vc Vc Vc] / 1E3;
S = S.' .* repmat( [dim dim], N, 1 );
nu = linspace(nu_0 - N * Ts, nu_f, N);

%% Results 
figure
hold on
stem(nu, dV_norm(1,:) * Vc, 'filled', 'r'); 
grid on;
ylabel('$\|\Delta \mathbf{V}\|_p$ [m/s]')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([nu(1) nu(end)])

siz = repmat(100, 1, 1);
figure 
view(3)
hold on
scatter3(S(1,7), S(1,8), S(1,9), siz, 'b', 'Marker', 'square');
scatter3(S(end,7), S(end,8), S(end,9), siz, 'b', 'Marker', 'o');
siz2 = repmat(100, sum(ti), 1);
scatter3( S(ti,7), S(ti,8), S(ti,9), siz2, 'Marker', 'x' );
plot3( S(:,7), S(:,8), S(:,9) ); 
legend('$\mathbf{s}(t_0)$', '$\mathbf{s}(t_f)$', '$\Delta \mathbf{V}_i$', '$\mathbf{s}(t)$', 'AutoUpdate', 'off');
hold off
grid on;
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
view(3)
hold on
scatter3( S(1,1) + S(1,7), S(1,2) + S(1,8), S(1,3) + S(1,9), siz, 'b', 'Marker', 'square' );
scatter3( S(end,1) + S(end,7), S(end,2) + S(end,8), S(end,3) + S(end,9), siz, 'b', 'Marker', 'o' );
scatter3( S(ti,1) + S(ti,1), S(ti,2) + S(ti,8), S(ti,3) + S(ti,9), siz2, 'Marker', 'x' );
plot3( S(:,1) + S(:,7), S(:,2) + S(:,8), S(:,3) + S(:,9) ); 
plot3( S(:,1), S(:,2), S(:,3) );
legend('$\mathbf{r}_c(t_0)$', '$\mathbf{r}_c(t_f)$', '$\Delta \mathbf{V}_i$', '$\mathbf{r}_c(t)$', '$\mathbf{r}_t(t)$', 'AutoUpdate', 'off');
hold off
grid on;
xlabel('$X$ [km]')
ylabel('$Y$ [km]')
zlabel('$Z$ [km]')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary function
% Newton equations of the co-orbital CR3BP
function [ds] = cr3bp_equations(mu, t, s, u)
    % Define the initial phase space vector
    r_t = s(1:3,:);                   % Synodic position vector
    x = s(1,:);                       % Synodic x coordinate
    y = s(2,:);                       % Synodic y coordinate 
    V = s(4:6,:);                     % Synodic velocity vector
    
    % Relevant system parameters
    mup(1) = 1 - mu;                                   % First primary normalized position
    mup(2) = mu;                                       % Second primary normalized position
    Rp(:,1) = [-mu; 0; 0];                             % Position vector of the first primary
    Rp(:,2) = [1-mu; 0; 0];                            % Position vector of the second primary
    
    r(1:3,:) = r_t(1:3,:) - Rp(:,1);                   % Synodic relative position of the target to the first primary
    r(4:6,:) = r_t(1:3,:) - Rp(:,2);                   % Synodic relative position of the target to the second primary

    R(1,:) = sqrt( dot(r(1:3,:), r(1:3,:), 1) );       % Distance to the first primary
    R(2,:) = sqrt( dot(r(4:6,:), r(4:6,:), 1) );       % Distance to the secondary primary
    
    % Compute the time flow of the system
    gamma = [x; y; zeros(1,size(x,2))];                % Inertial acceleration terms
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds = [V; gamma]; 

    % Gravitational forces
    Accg = - mup(1) ./ R(1,:).^3 .* r(1:3,:) - mup(2) ./ R(2,:).^3 .* r(4:6,:);
    ds(4:6,:) = ds(4:6,:) + Accg;
    
    % Compute the time flow of the system
    r_r = s(7:9,:);                                    % Synodic position vector
    x = s(7,:);                                        % Synodic x coordinate
    y = s(8,:);                                        % Synodic y coordinate 
    V = s(10:12,:);                                    % Synodic velocity vector

    gamma = [x; y; zeros(1,size(x,2))];                % Inertial acceleration terms
    gamma = gamma + [0 2 0; -2 0 0; 0 0 0] * V;
    ds(7:12,:) = [V; gamma]; 
 
    % Gravitational forces
    F =   + mup(1) * ( r(1:3,:) ./ R(1,:).^3 - (r_r + r(1:3,:)) ./ sqrt( dot(r_r + r(1:3,:), r_r + r(1:3,:), 1) ).^3 );
    F = F + mup(2) * ( r(4:6,:) ./ R(2,:).^3 - (r_r + r(4:6,:)) ./ sqrt( dot(r_r + r(4:6,:), r_r + r(4:6,:), 1) ).^3 );
 
    % Control force 
    ds(10:12,:) = ds(10:12,:) + F;
end

% Compute the STM of the relative CRB3P
function [Phi] = wrapper_STM(c2, nu_0, nu_f, N)
    % Pre-allocation
    Phi = zeros(6, 6 * N);
    STM = zeros(6, 6 * N);

    nu = linspace(nu_0, nu_f, N);
    
    for i = 1:length(nu)
        idx = 1 + 6 * (i-1) : 6 * i;
        dt = nu(i) - nu(1);
        STM(:,idx) = CR3BP_STM(c2, dt);
        Phi(:,idx) = STM(:,idx);
    end
end

% STM of the relative CR3BP
function [Phi] = CR3BP_STM(c2, delta_t)
    % State space matrix 
    Omega = 2 * [0 1 0; -1 0 0; 0 0 0];             % Coriolis term
    H = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];            % Hessian of the Hamiltonian
    A = [zeros(3) eye(3); H Omega];                 % State space matrix
    
    Phi = expm(A * delta_t);
end