%% Optimal Linear Rendezvous via ADMM in the CR3BP %% 
% Sergio Cuevas del Valle
% Date: 17/09/25
% File: Serra2018.m 
% Issue: 0 
% Validated: 

%% C3BP rendezvous, Serra 2018 %% 
% Solve for the time-fixed optimal rendezvous in the CR3BP (EM L2) via ADMM
% and PVT %

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Parameters
Lc = 384399e3;          % Characteristic length
Tc = 2.361e6;           % Characteristic time
mu = 0.0121505856;      % Gravitational constan of the system

% Dimensionalization (canonical units)
Vc = Lc / Tc * 2*pi;    % Characteristic velocity 

% Mission time
nu_0 = 3.322;           % Initial clock
nu_f = 4.737;           % Final clock

% Initial relative conditions 
x0 = [6449.40 65117.03 22814.91 -0.0312 0.0392 0.2114];    % Initial conditions
xf = [59066.09 67728.64 84015.47 -0.1087 0.1616 -0.1730];  % Final conditions

x0 = x0 ./ [Lc Lc Lc Vc Vc Vc];
xf = xf ./ [Lc Lc Lc Vc Vc Vc];

x0 = x0.'; 
xf = xf.';

% Number of possible impulses 
N = 30;

%% Define the rendezvous problem and the STM %%
% Time span
nu = linspace(0, 1, N);

% Control input matrix
B = [zeros(3); eye(3)];
B = repmat( B, 1, length(nu) );

c2 = 3.190425213622208;

% Pre-allocation
STM = @(nu_f)wrapper_STM(c2, nu_0, nu_f, N);

%% Final mission definition 
K = Inf;                                                       % Maximum number of impulses
myMission = Missions.TimeMission(nu, nu_0, STM, B, x0, xf, K); % Mission

%% Thruster definition 
dVmin = 0;                                                     % Minimum control authority
dVmax = 0.8;                                                   % Maximum control authority
myThruster = Actuator(src.VectorNorm.L2, dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myDualProblem   = MinimumTimeSolvers.NeustadtSolver(myMission, myThruster);

iter = 1;                           % Number of interations
time = zeros(1,iter);               % Computational cost
dV = zeros(3 * 1, N);               % Impulses of the two algorithms

% Optimization
rho = 1/N;                          % AL parameter 
eps = [1e-6; 1e-5];                 % Numerical tolerance

for i = 1:iter
    % Dual resolution
    [~, sol, ~, tf, myDualProblemSolved] = myDualProblem.Solve(eps, rho^(3/2));
    time(1,i) = myDualProblemSolved.SolveTime;

    lambda = reshape(sol(1:6), 6, []);
    p = reshape(sol(7:end), 3, []);
    dV(1:3,:) = myDualProblemSolved.u;
end

%% Outcome
% Cost function
dV_norm(1,:) = myThruster.p.ComputeVectorNorm( dV(1:3,:) );

% Norm of the primer vector 
p_norm = myThruster.q.ComputeVectorNorm( p );

% Impulsive times 
ti(1,:) = dV_norm(1,:) ~= 0;

% Results
cost(1) = sum(dV_norm(1,ti(1,:)), 2) * Vc;

Nopt = sum(ti, 2);
Time = mean(time, 2);
error(1,:) = sqrt( dot(myPrimalProblemSolved.e, myPrimalProblemSolved.e, 1) );

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(nu), 6 * 1);
s(1,:) = [x0.'];

% Computation
for i = 1:length(nu)
    for j = 1:1
        % Propagate 
        if (i > 1)
            prev_idx = 1 + 6 * (i - 2) : 6 * (i - 1);
            curr_idx = 1 + 6 * (i - 1) : 6 * (i - 0);
    
            Phi1 = reshape(STM(:,prev_idx), [6 6]);
            Phi2 = reshape(STM(:,curr_idx), [6 6]);
    
            state_idx = 1 + 6 * (j-1) : 6 * j;
            s(i,state_idx) = s(i-1,state_idx) * (Phi2 * Phi1^(-1)).';
        end
    
        % Add maneuver
        cntrl_indx = 4 + 6 * (j - 1) : 4 + 6 * (j - 1) + 2;
        plan_idx = 1 + 3 * (j-1) : 3 * j;
        s(i,cntrl_indx) = s(i,cntrl_indx) + dV(plan_idx,i).';
    end
end

% Dimensionalization 
dim = [Lc Lc Lc Vc Vc Vc];
s = s .* repmat(dim, N, 1) / 1e3;

%% Save results 
save +Papers_EuroGNC_2026\+RendezvousCR3BP\MinTimeResultsSerraL2

%% Results 
% Norm of the primer vector
nu_imp = nu( ti(1,:) );

figure
hold on
scatter(nu_imp, ones(1, length(nu_imp)), 1e2, 'r', 'Marker', 'x')
legend('$t_i$', 'AutoUpdate', 'off')
plot(nu, p_norm, 'b');
yline(1, '--')
grid on;
ylabel('$\|\mathbf{p}\|_q$')
xlabel('$t$')
xlim([nu(1) nu(end)])

figure
hold on
stem(nu, dV_norm(1,:) * Vc, 'filled', 'r'); 
legend('Neustadt', 'AutoUpdate', 'off')
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
scatter3(s(1,1), s(1,2), s(1,3), siz, 'b', 'Marker', 'square');
scatter3(s(end,1), s(end,2), s(end,3), siz, 'b', 'Marker', 'o');


for j = 1:1
    % Plot each trajectory and the corresponding control law
    siz2 = repmat(100, sum(ti(j,:)), 1);
    state_idx = [1 2 3] + 6 * (j - 1);
    impulses = ti(j,:);

    scatter3( s(impulses,state_idx(1)), s(impulses,state_idx(2)), s(impulses,state_idx(3)), siz2, 'Marker', 'x' );
    plot3( s(:,state_idx(1)), s(:,state_idx(2)), s(:,state_idx(3)) ); 
end 
legend('$\mathbf{s}_0$', '$\mathbf{s}_f$', '$\Delta \mathbf{V}_i^N$', 'Neustadt', 'AutoUpdate', 'off');

hold off
grid on;
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));


%% Auxiliary function 
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

function [Phi] = CR3BP_STM(c2, delta_t)
    % State space matrix 
    Omega = 2 * [0 1 0; -1 0 0; 0 0 0];             % Coriolis term
    H = [1+2*c2 0 0; 0 1-c2 0; 0 0 -c2];            % Hessian of the Hamiltonian
    A = [zeros(3) eye(3); H Omega];                 % State space matrix
    
    Phi = expm(A * delta_t);
end