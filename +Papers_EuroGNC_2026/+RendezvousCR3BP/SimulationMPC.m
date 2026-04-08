%% Optimal Linear Rendezvous via ADMM in the CR3BP %% 
% Sergio Cuevas del Valle
% Date: 14/10/25
% File: SimulationMPC.m 
% Issue: 0 
% Validated: 

%% CR3BP rendezvous via MPC % 
% Solve for a time-fixed optimal rendezvous in the CR3BP (EM L2) via ADMM, MPC and PVT %

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
Lc = 384399E3;          % Characteristic length [m]
Tc = 2.361E6;           % Characteristic time [s]
mu = 0.0121505856;      % Gravitational constan of the system
c2 = 3.190425213622208; % Richardson coefficient (L1)
c2 = 4.284390615378499; % Richardson coefficient (L2)
    
% Dimensionalization (canonical units)
Vc = Lc / Tc * 2*pi;    % Characteristic velocity 

% Mission time
TargetPeriod = 2.761104629643622; 
ChaserPeriod = 2.754937512093445; 

nu_0 = 0;               % Initial clock
nu_f = ChaserPeriod;    % Final clock

% Initial target conditions 
x0 = [0.824024728136525; 0; -0.054501847320725; 0; 0.164671964079122; 0];          

% Initial chaser conditions
xc = [0.823639438925721; 0; +0.043281569720089; 0; 0.152567980892620; 0];

% Target orbital propagation 
options = odeset('RelTol', 2.25E-14, 'AbsTol', 1E-22);      % Numerical integration
tspan = [nu_0 nu_f];
[~, St] = ode113(@(t,s)cr3bp_equations(mu, t, s, zeros(3,1)), tspan, [x0; zeros(6,1)], options);

% Chaser orbital propagation 
[~, Sc_free] = ode113(@(t,s)cr3bp_equations(mu, t, s, zeros(3,1)), tspan, [xc; zeros(6,1)], options);

% Co-orbital initial conditions
x0 = [x0; xc - x0];

% Number of possible impulses 
N  = 100;                       % Number of steps
Ts = (nu_f - nu_0) / N;         % Sampling time

% Model to be used
model = 2;                      % Use RLLM model (1) or account for the target's trajectory (any other value)

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
% Optimization
eps = 1E-5;                % Numerical tolerance
dV_final = zeros(3,N);     % Maneuver sequence
S = zeros(12,N);           % Realized trajectory
OptTime = zeros(1, N);     % Optimization time  
iter = 1;                  % Iteration index
Ninit = N;                 % Initial number of problems
m = 6;                     % Dimension of the co-orbital state

% Pre-allocation of the impulses
dV = zeros(3,N);       
initial_guess = [];

while ( iter <= Ninit )
    if ( N > 1 )
        % Update of the mission
        nu = linspace(nu_0, nu_f, N);
        
        % STM
        if ( model == 1 )
            % Integrate the RLLM variational model
            Phi = wrapper_STM(c2, nu, N);

        else
            % Integrate the RLM variational model
            Phi = reshape(eye(6), [], 1);
            [~, s] = ode113( @(t,s)cr3bp_var( mu, t, s, zeros(3,1) ), nu, [x0(1:6); Phi], options );

            if ( N == 2 )
                Phi = [eye(6) reshape( s(end,7:end), 6, [] )];
            else
                Phi = reshape( s(:,7:end).', 6, [] );
            end
        end

        Binp = repmat( B, 1, N );
        myMission = Missions.FuelMission(nu, Phi, Binp, x0(7:end,1), zeros(6,1), K);
           
        % Define the ADMM problem 
        rho = 1 / N;                       % AL parameter 
%         myDualProblem = MinimumNormSolvers.NeustadtSolver(myMission, myThruster);
        myPrimalProblem = MinimumNormSolvers.PrimalSolver(myMission, myThruster);

        % Optimization
%         [~, sol, ~, myDualProblemSolved] = myDualProblem.Solve(eps, rho^(3/2), 1);
        [~, sol, ~, myPrimalProblemSolved] = myPrimalProblem.Solve( 1/rho );

        OptTime(iter) = myPrimalProblemSolved.SolveTime;
    
        if ( ~isempty( sol ) )
            % New maneuver sequence
            dV = myPrimalProblemSolved.u;

            % New initial guess
%             initial_guess.x = sol([1:m m+4:end]);
%             initial_guess.z = initial_guess.x;
        else
            % New maneuver sequence
            dV = dV(:,2:end);

            % New initial guess
            if ( ~isempty(initial_guess) )
%                 initial_guess.x = initial_guess.x([1:m m+4:end]);
%                 initial_guess.z = initial_guess.x;
            end
        end
    else
        dV = dV(:,end);
    end
    
    % Apply the impulse 
    x0(10:12,1)      = x0(10:12,1) + dV(:,1);
    S(:,iter)        = x0;
    dV_final(:,iter) = dV(:,1);

    % Integration the coasting solution
    tspan = [nu_0 nu_0 + Ts];
    [~, s] = ode113( @(t,s)cr3bp_equations(mu, t, s, zeros(3,1)), tspan, x0, options );
    
    % Update the IVP
    nu_0 = nu_0 + Ts;           % Update the initial clock
    x0 = s(end,:).';            % New relative initial conditions

    % Update of the MPC loop 
    N = N - 1;                  % Reduce the number of optimization steps

    fprintf("Iteration: %d. Converged? %d\n", iter, ~isempty( sol ))
    iter = iter + 1;
end

%% Outcome
% Cost function
dV_norm = myThruster.p.ComputeVectorNorm( dV_final(1:3,:) ); 
ti = dV_norm ~= 0;
cost = sum( dV_norm, 2 ) * Vc;
Nopt = sum(ti,2);

%% Save results 
% save +Papers_EuroGNC_2026\+RendezvousCR3BP\MPC_L2_N100

%% Dimensionalizations 
dim = [Lc Lc Lc Vc Vc Vc] / 1E3;

St     = St        .* repmat( [dim dim], size(St,1), 1 );
Sc_free = Sc_free  .* repmat( [dim dim], size(Sc_free,1), 1 );

N = size(S,2);
S  = S.' .* repmat( [dim dim], N, 1 );

% Absolute trajectories 
Sc = S(:,1:6) + S(:,7:12);

% Timing 
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

figure
hold on
plot(OptTime, 'ro-'); 
grid on;
ylabel('Opt. time [s]')
xlabel('$N$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
xlim([1 Ninit])

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
%%
siz = repmat(100, 1, 1);
siz2 = repmat(100, sum(ti), 1);
figure 
view(3)
hold on
plot3( St(:,1), St(:,2), St(:,3) );
plot3( Sc_free(:,1), Sc_free(:,2), Sc_free(:,3) ); 
scatter3( Sc(1,1), Sc(1,2), Sc(1,3), siz, 'b', 'Marker', 'square' );
scatter3( Sc(end,1), Sc(end,2), Sc(end,3), siz, 'b', 'Marker', 'o' );
plot3( Sc(:,1), Sc(:,2), Sc(:,3) ); 
legend('$\mathbf{r}_t(t)$', '$\mathbf{r}_c(t)$', '$\mathbf{r}_c(t_0)$', '$\mathbf{r}_c(t_f)$', '$\mathbf{r}_c^u(t)$', 'AutoUpdate', 'off');
hold off
grid on;
xlabel('$X$ [km]')
ylabel('$Y$ [km]')
zlabel('$Z$ [km]')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

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
    ds(10:12,:) = ds(10:12,:) + F + u;
end

% Rereference target dynamics + variational equations
function [ds] = cr3bp_var(mu, t, s, u)
    % Define the initial phase space vector
    r_t = s(1:3,:);                                    % Synodic position vector
    x = s(1,:);                                        % Synodic x coordinate
    y = s(2,:);                                        % Synodic y coordinate 
    V = s(4:6,:);                                      % Synodic velocity vector
    
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
    ds(4:6,:) = ds(4:6,:) + Accg + u;
    
    % Define the Jacobian of the co-orbital model 
    Sigma = [0 2 0; -2 0 0; 0 0 0];
    eps = [r(1:3,:) ./ R(1,:) r(4:6,:) ./ R(2,:)];
    kappa = [mup(1) ./ R(1,:).^3 mup(2) ./ R(2,:).^3];
    H = -sum(kappa) * eye(3) + 3 * kappa(1) * eps(:,1) * eps(:,1).' + 3 * kappa(2) * eps(:,2) * eps(:,2).';
    J = [zeros(3) eye(3); H Sigma];

    % Differential system 
    Phi = reshape(s(7:end), 6, 6);
    dJ = J * Phi; 
    ds = [ds; reshape(dJ, [], 1)];
end

% Compute the STM of the relative CRB3P
function [Phi] = wrapper_STM(c2, nu, N)
    % Pre-allocation
    Phi = zeros(6, 6 * N);
    STM = zeros(6, 6 * N);
    
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