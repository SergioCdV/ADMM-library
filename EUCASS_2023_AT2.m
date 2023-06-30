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
I = diag([1 1 1]);   % Inertia matrix

%% Define the rendezvous problem and the STM 
tf = 600;             % Time of flight [s]
    
% Time span
t = linspace(0, tf, 100);
dt = diff(t); 
dt = dt(1);

% Relative initial conditions
x0 = [zeros(3,1); 1; zeros(3,1)];  

%% Optimization
% Define the ADMM problem 
rho = 1;        % AL parameter 

% Integration options
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25e-14);
    
% Preallocation 
time = zeros(1,length(t));  % Computational time
U = zeros(3,length(t));     % Control law
S = zeros(length(t), 7);    % State trajectory
S(1,:) = x0.';              % Initial conditions
Bc = [zeros(3); I^(-1)];    % Control input matrix

% Reference state
ref = [0; 0; sind(45/2); cosd(45/2); zeros(3,1)];  

[~, saux] = ode113(@(t,s)dynamics(I, t, s, zeros(3,1)), t, S(1,:), options);
sf = saux(end,:);

% MPC
for i = 1:length(t)-1
    % Optimization & planning
    n = 3 * length(t(i:end));
    m = n;
    A = eye(n);
    B = -eye(n);        
    c = zeros(m,1);
    
    p = repmat(3, 1, n/3);
    cum_part = cumsum(p);

    Obj = @(x,z)(objective(cum_part, z));
    X_update = @(x,z,u)(x_update(dt, t(i:end), S(i,:).', ref, I, Bc, x, z, u));
    Z_update = @(x,z,u)(z_update(rho, p, x, z, u));
        
    Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
    Problem.alpha = 1;
    Problem.QUIET = false; 

    tic
    [x, z, Output] = Problem.solver();
    T = reshape(x(:,end), 3, []);
    time(i) = toc; 

    % Save the slew
    U(:,i) = T(:,1);

    % Integration of the dynamics
    [~, saux] = ode113(@(t,s)dynamics(I, t, s, U(:,i)), [0 dt], S(i,:), options);
    S(i+1,:) = saux(end,:);
end

%% Outcome
R = zeros(size(S,1), 4);
for i = 1:size(S,1)
    R(i,1:4) = S(i,1:4);
end

cost_admm = sum( sqrt(dot(T, T, 1)) );

%% Results 
figure
hold on
stem(t, sqrt(dot(U,U,1)), 'filled'); 
grid on;
ylabel('$\tau$')
xlabel('$t$')
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));

plotTripleEvolution(t, S);

figure
view(3)
hold on
no = [eye(3); zeros(1,3)];
quiver3(zeros(3,1), zeros(3,1), zeros(3,1), no(1:3,1), no(1:3,2), no(1:3,3), 'b', 'LineWidth', 1.0); 

nref = no;
for i = 1:1
    for j = 1:size(no,2)
        aux = quaternion_product(no(:,j), quaternion_inverse(ref(1:4)));
        nref(:,size(no,2)*(i-1)+j) = quaternion_product(ref(1:4), aux);
    end
end
quiver3(zeros(3,1), zeros(3,1), zeros(3,1), nref(1:3,1), nref(1:3,2), nref(1:3,3), 'r', 'LineWidth', 1.0);

aux = quaternion_product([0 1 0 0].', quaternion_inverse(R(1,1:4)));
C(1,:) = quaternion_product(R(1,1:4).', aux).';
for i = 2:size(S,1)
    aux = quaternion_product([0 1 0 0].', quaternion_inverse(R(i,1:4)));
    C(i,:) = quaternion_product(R(i,1:4).', aux).';
end
plot3(C(:,1), C(:,2), C(:,3), 'k');
grid on;

legend('$\mathcal{P}_0$', '$\mathcal{P}_{ref}$', '$\mathcal{C}$', 'AutoUpdate', 'off', 'Location', 'best')

m = 100;
[aa, bb, cc] = sphere(m);
h = surf(aa, bb, cc);
set(h, 'FaceAlpha', 0.1)
shading interp

hold off
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary functions 
% Dynamics 
function [ds] = dynamics(I, t, s, tau)
    % Initial variables 
    q = s(1:4);             % Error quaternion 
    omega = s(5:7);         % Error angular velocity 

    % Kinematics 
    dq = 0.5 * quaternion_product([omega;0], q);

    % Dynamics 
    domega = I^(-1) * (tau - cross(omega, I * omega));

    % Final vector field 
    ds = [dq; domega];
end

% Hat map
function [S] = hat_map(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% Quaternion matrix 
function [Q] = quaternion_matrix(dq)
    S = hat_map(dq(1:3)); 
    Q = [dq(4) * eye(3)-S dq(1:3); -dq(1:3).' dq(4)];
end

% Quaternion product 
function [q2] = quaternion_product(dq, q1)
    Q = quaternion_matrix(dq);
    q2 = Q * q1;
end

% Quaternion inverse 
function [q2] = quaternion_inverse(q1)
    q2(1:3,1) = -q1(1:3);
    q2(4,1) = q1(4);
end

% Exponential mapping 
function [dq] = exp_map(v)
    if (norm(v) ~= 0)
        dq = [v/norm(v) * sin(norm(v)); cos(norm(v))];
    else
        dq = [0;0;0;1];
    end
end

% Proximal minimization for X
function [x] = x_update(dt, t, s, ref, I, B, x, z, u) 
   % Computation of the STM 
   STM = zeros(length(t), 36);
   STM(1,:) = reshape(eye(6), 1, []);

   G = zeros(6, 3 * length(t));
   omega = s(5:7);

   for i = 2:length(t)
       % STM
       A = [-hat_map(omega) 0.5 * eye(3); zeros(3) I^(-1) * hat_map(I * omega)];
       Omega = expm(A * dt * (i-1));
       STM(i,:) = reshape(Omega, 1, []);
   end

  M = reshape(STM(end,:), [6 6]);
  for i = 1:length(t)
       Phi = M * reshape(STM(i,:), [6 6])^(-1);
       AB = dt/2 * (eye(6) + Phi) * B;
       G(:,1+3*(i-1):3*i) = AB;
  end

   % Computation of the missvector 
   sf = quaternion_product(exp_map(omega * (t(end)-t(1))/2), s(1:4,1));
   e(1:4,1) = quaternion_product(ref(1:4), quaternion_inverse(sf(1:4,1)));
   e = e(1:3);
   e(4:6) = ref(5:7)-omega; 

   % Sequence update
   pInvA = (eye(size(G,2))-pinv(G)*G);
   Atb = pinv(G) * e;
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
    ylabel('$\mathbf{q}$');
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
    ylabel('$\mathbf{\omega}$');
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