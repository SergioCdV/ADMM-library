%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: EUCASS_2023_AT2.m 
% Issue: 0 
% Validated: 

%% Regulation / attitude slew %% 
% Solve for the time-fixed attitude slew optimal L2 problem using ADMM 

close; 
clear; 
clc

set_graphics();

%% Define the target orbit
% Parameters
I = diag([1 2 3]);   % Inertia matrix

%% Define the rendezvous problem and the STM 
tf = 60;             % Time of flight [s]
    
% Time span
t = linspace(0, tf, 500);
dt = diff(t); 
dt = dt(1);

% Relative initial conditions
x0 = [sqrt(2)/2; 0; 0; sqrt(2)/2; zeros(3,1)];  

% Relative final conditions 
xf = [zeros(3,1); 1; zeros(3,1)];

%% Optimization
% Define the ADMM problem 
rho = 5e2;        % AL parameter 

% Constraint 
[m,n] = size(Phi); 
A = eye(n);
B = -eye(n);        
c = zeros(m,1);

p = repmat(3,1,length(t));
cum_part = cumsum(p);

% Integration options
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25e-14);
    
% Preallocation 
time = zeros(1,length(t));  % Computational time
U = zeros(3,length(t));     % Control law
S = zeros(length(t), 7);    % State trajectory
S(1,:) = x0.';              % Initial conditions

% MPC
for i = 1:length(t)-1
    % Optimization & planning
    Obj = @(x,z)(objective(cum_part, z));
    X_update = @(x,z,u)(x_update(pInvA, Atb, x, z, u));
    Z_update = @(x,z,u)(z_update(rho, p, x, z, u));
        
    Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
    Problem.alpha = 1;
   
    tic
    [x, z, Output] = Problem.solver();
    T = reshape(x(:,end), 3, []);
    time(i) = toc; 

    % Integration of the dynamics
    [~, saux] = ode113(@(t,s)dynamics(I, t, s, T(1,:)) , [0 1e-3], S(i,:), options);
    [~, saux] = ode113(@(t,s)dynamics(I, t, s, zeros(3,1)), [0 dt], saux(end,:), options);

    S(i+1,:) = saux(end,:);

    % Save the slew
    U(i,:) = T(1,:);
end

cost_admm = sum( sqrt(dot(T, T, 1)) );

%% Results 
figure
hold on
stem(t, sqrt(dot(T,T,1)), 'filled'); 
grid on;
ylabel('$\Delta V$')
xlabel('$t$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot3(s(:,1), s(:,2), s(:,3)); 
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

plotTripleEvolution(t, S);

%% Auxiliary functions 
% Dynamics 
function [ds] = dynamics(I, t, s, alpha)
    % Initial variables 
    q = s(1:4);             % Error quaternion 
    omega = s(5:7);         % Error angular velocity 

    % Kinematics 
    dq = 0.5 * quaternion_product([omega;0], q);

    % Dynamics 
    domega = alpha - I^(-1) * cross(omega, I * omega);

    % Final vector field 
    ds = [dq; domega];
end

% Quaternion product 
function [q2] = quaternion_product(dq, q1)
    S = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0]; 
    Q = [dq(4) * eye(3)-S dq(1:3); -dq(1:3).' dq(4)];
    q2 = Q * q1;
end

% Proximal minimization for X
function [x] = x_update(pInvA, Atb, x, z, u) 
   x = pInvA*(z-u)+Atb;
end

% Proximal minimization for Z
function [z] = z_update(rho, p, x, z, u)
    % cumulative partition
    cum_part = cumsum(p);

    % z-update
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:cum_part(i);
        z(sel) = shrinkage(x(sel) + u(sel), 1/rho);
        start_ind = cum_part(i) + 1;
    end
end

% Proximal minimization of the l2 norm
function z = shrinkage(x, kappa)
    z = max(0, 1 - kappa/norm(x))*x;
end

% Objective function
function [p] = objective(cum_part, z)
    obj = 0;
    start_ind = 1;
    for i = 1:length(cum_part)
        sel = start_ind:cum_part(i);
        obj = obj + norm(z(sel));
        start_ind = cum_part(i) + 1;
    end
    p = ( obj );
end

%Function to plot the three relative state evolutions in the same figure 
function plotTripleEvolution(tspan_misg, St_misg)    
    figure
    subplot(1,2,1)
    hold on
    for i = 1:4
        % Configuration space evolution
        plot(tspan_misg, St_misg(:,i));
    end
    legend('$q_x$', '$q_y$', '$q_z$', '$q_w$', 'AutoUpdate', 'off');
    hold off;
    xlabel('$t$');
    ylabel('$\mathbf{\rho}$');
    grid on;

    subplot(1,2,2)
    hold on
    for i = 1:3
        % Configuration space evolution
        plot(tspan_misg, St_misg(:,4+i));
    end
    legend('$\omega_x$', '$\omega_y$', '$\omega_z$', 'AutoUpdate', 'off');
    hold off;
    xlabel('$t$');
    ylabel('$\mathbf{\omega}}$');
    grid on;
end

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