%% Optimal Linear Slew via ADMM %% 
% Sergio Cuevas del Valle
% Date: 06/10/25
% File: AxisymmetricSlew.m 
% Issue: 0 
% Validated: 

%% Axisymmetric attitude slew %% 
% Solve for a minimum time attitude slew via prox optimization %

close; 
clear; 
clc

set_graphics();

%% Define the target orbit and mission parameters
% Mission time
nu_0 = 0;               % Initial clock

% Initial conditions 
x0 = [0 0 0 1 1 0 0];                  % Initial conditions (attitude + angular velocity)
xf = [sqrt(2)/2 0 0 sqrt(2)/2 0 0 0];  % Final conditions   (attitude + angular velocity)

% Problem parameters 
I = diag( [1 1 3] );                   % Inertia tensor of the spacecraft

b = I(3,3) / I(1,1) - 1;               % Oblateness parameter    
omega0 = norm( x0(5:7) );              % Initial angular velocity
r0 = x0(end);                          % Initial angular velocity along the last axis
    
% Number of possible impulses 
N = 100;

% Compute the STM across a given time span 
nu_f = 100; 
nu = linspace(nu_0, nu_f, N);

STM = zeros(7, 7 * N);

for i = 1:N 
    idx = 1 + 7 * (i-1) : 7 * i;
    delta_t = nu(i) - nu(1);
    STM(:,idx) = Paper_Attitude_EuroGNC_2025.SlewSTM(b, omega0, r0, delta_t);
end

%% Define the rendezvous problem and the STM %%
% Control input matrix
B = [zeros(4); eye(3)];
B = repmat( B, 1, length(nu) );

%% Final mission definition 
K = Inf;                                                % Maximum number of impulses
myMission = LinearMission(nu, Phi, B, x0, xf, K);       % Mission

%% Thruster definition 
dVmin = 0;                                              % Minimum control authority
dVmax = Inf;                                            % Maximum control authority
myThruster = thruster('L1', dVmin, dVmax);

%% Optimization
% Define the ADMM problem 
myDualProblem = MinimumTimeSolvers.NeustadtSolver(myMission, myActuator);

iter = 25;                          % Number of interations
time = zeros(2,iter);               % Computational cost
dV = zeros(3 * 2, N);               % Impulses of the two algorithms

% Optimization
rho = 1/N;                          % AL parameter 
eps = [1e-6; 1e-5];                 % Numerical tolerance

for i = 1:iter
    % Dual resolution
    [~, sol, ~, myDualProblemSolved] = myDualProblem.Solve(eps, rho^(3/2));
    time(1,i) = myDualProblemSolved.SolveTime;

    lambda = reshape(sol(1:6), 6, []);
    p = reshape(sol(7:end), 3, []);
    dV(1:3,:) = myDualProblemSolved.u;
end

%% Outcome
% Cost function
switch (myActuator.p)
    case 'L1'
        dV_norm(1,:) = sum( abs(dV(1:3,:) ), 1);

    case 'L2'
        dV_norm(1,:) = sqrt( dot(dV(1:3,:), dV(1:3,:), 1) );

    case 'Linfty'
        dV_norm(1,:) = max( abs( dV(1:3,:) ) );
end

% Norm of the primer vector 
switch (myActuator.q)
    case 'L1'
        p_norm = sum( abs(p), 1 );
    case 'L2'
        p_norm = sqrt( dot(p, p, 1) );
    case 'Linfty'
        p_norm = max( abs(p) );
end

% Impulsive times 
ti(1,:) = dV_norm(1,:) ~= 0;

% Results
cost(1) = sum(dV_norm(1,ti(1,:)), 2);
Nopt = sum(ti, 2);
Time = mean(time, 2);
error(1,:) = sqrt( dot(myDualProblemSolved.e, myDualProblemSolved.e, 1) );

%% Chaser orbit reconstruction 
% Preallocation 
s = zeros(length(nu), 6);
s(1,:) = [x0.'];

% Computation
for i = 1:length(nu)
    % Propagate 
    if (i > 1)
        prev_idx = 1 + 6 * (i - 2) : 6 * (i - 1);
        curr_idx = 1 + 6 * (i - 1) : 6 * (i - 0);

        Phi1 = reshape(STM(:,prev_idx), [6 6]);
        Phi2 = reshape(STM(:,curr_idx), [6 6]);

        state_idx = 1:6;
        s(i,state_idx) = s(i-1,state_idx) * (Phi2 * Phi1^(-1)).';
    end
    
    % Add maneuver
    cntrl_indx = 4:6;
    plan_idx = 1:3;
    s(i,cntrl_indx) = s(i,cntrl_indx) + dV(plan_idx,i).';
end

%% Save results 
save +Paper_Attitude_2025\MinTimeL2

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
stem(nu, dV_norm(1,:), 'filled', 'r'); 
grid on;
ylabel('$\|\mathbf{u}\|_p$ [m/s]')
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

% Plot each trajectory and the corresponding control law
siz2 = repmat(100, sum(ti(1,:)), 1);
state_idx = [1 2 3];
impulses = ti(1,:);

scatter3( s(impulses,state_idx(1)), s(impulses,state_idx(2)), s(impulses,state_idx(3)), siz2, 'Marker', 'x' );
plot3( s(:,state_idx(1)), s(:,state_idx(2)), s(:,state_idx(3)) ); 

legend('$\mathbf{s}_0$', '$\mathbf{s}_f$', '$\mathbf{u}_i$', 'AutoUpdate', 'off');

hold off
grid on;
xlabel('$x$ [-]')
ylabel('$y$ [-]')
zlabel('$z$ [-]')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));
