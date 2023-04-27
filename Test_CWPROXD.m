%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% Test exercise VIII %% 
% Solve for the time-fixed CW optimal L2 problem using ADMM and PVT

close; 
clear; 
clc

%% Define the rendezvous problem and the STM 
% Parametric study 
parametric = false; 

if (parametric)
    TOF = linspace(0.1,50,100); 
else
    TOF = pi;
end

ratio = zeros(1,length(TOF));
cost = zeros(1,length(TOF));

%% Optimization
for j = 1:length(TOF)
    % Time span
    t = linspace(0,TOF(j),100);
    n = 1;
    
    % Relative initial conditions 
    x0 = [1; 2; 3; -1; -pi; 2]; 
    
    % Relative final conditions 
    xf = zeros(6,1);
    
    % CW STM
    STM(:,1:6) = [4-3*cos(n*t.') -6*n*t.'+6*sin(n*t.') zeros(length(t),1) 3*n*sin(n*t.') -6*n+6*n*cos(n*t.') zeros(length(t),1)]; 
    STM(:,7:12) = [zeros(length(t),1) ones(length(t),1) zeros(length(t),4)];
    STM(:,13:18) = [zeros(length(t),2) cos(n*t.') zeros(length(t),2) -n*sin(n*t.')];
    STM(:,19:24) = [sin(n*t.')/n -2/n+2*cos(n*t.')/n zeros(length(t),1) cos(n*t.') -2*sin(n*t.') zeros(length(t),1)];
    STM(:,25:30) = [2/n-2*cos(n*t.')/n 4/n*sin(n*t.')-3*t.' zeros(length(t),1) 2*sin(n*t.') -3+4*cos(n*t.') zeros(length(t),1)];
    STM(:,31:36) = [zeros(length(t),2) sin(n*t.')/n zeros(length(t),2) cos(n*t.')];
    
    Phi = zeros(6, 3 * length(t));
    M = reshape(STM(end,:), [6 6]);
    for i = 1:length(t)
        Phi(:,1+3*(i-1):3*i) = M * reshape(STM(i,:), [6 6])^(-1) * [zeros(3); eye(3)];
    end
    
    % Residual or missvector
    b = xf-M*x0;

    % Maximum L2 norm 
    dVmax = 20;
    
    %% Define the ADMM problem 
    rho = 10;        % AL parameter 
    
    % pre-factor
    Atb = pinv(Phi)*b;
    pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);
    
    p = repmat(3,1,length(t));
    cum_part = cumsum(p);
        
    % Constraint 
    [m,n] = size(Phi); 
    A = eye(2 * n);
    B = -eye(2 * n);        
    c = zeros(m,1);

    % Create the functions to be solved 
    Obj = @(x,z)(objective(cum_part, z));
    X_update = @(x,z,u)(x_update(pInvA, Atb, n, x, z, u));
    Z_update = @(x,z,u)(z_update(Phi, b, dVmax, cum_part, n, rho, x, z, u));
    
    % Problem
    Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
    Problem.alpha = 1;
    
    %% Optimization
    [x, z, Output] = Problem.solver();
    dV = reshape(x(1:n,end), 3, []);
    p = reshape(x(n+1:2*n,end), 3, []);

    cost_admm = sum(sqrt(dot(dV,dV,1)));

    %% Pruning
    [dV2, cost(j)] = PVT_pruner(STM, [zeros(3); eye(3)], dV, 'L2');

    %% Optimality (PVT) checks 
    G = zeros(6,6);
    for i = 1:length(t)
        D = M * reshape(STM(i,:), [6 6])^(-1) * [zeros(3); eye(3)];
        G = G + norm(dV(:,i)) * (D*D.');
    end

    lambda = G\b;

    mag_p = sqrt(dot(p,p,1));
    error = abs(1-mag_p);
   
    %% Outcome 
    res(:,1) = b-Phi*x(1:n,end);
    res(:,2) = b-Phi*reshape(dV2, [], 1);
    ratio(j) = 1-cost(j)/cost_admm;
end

%% Results 
figure
hold on
plot(1:Output.Iterations, Output.objval); 
yline(cost, '--')
legend('ADMM cost', 'PVT cost')
grid on;
ylabel('$\Delta V_T$')
xlabel('Iteration $i$')

figure
hold on
yline(dVmax, '--');
stem(t, sqrt(dot(dV,dV,1)), 'filled'); 
grid on;
legend('$\Delta V_{max}$', '$\Delta V_{ADMM}$')
ylabel('$\Delta V$')
xlabel('$t$')

figure
hold on
yline(dVmax, '--');
stem(t, sqrt(dot(dV2,dV2,1)), 'filled'); 
grid on;
legend('$\Delta V_{max}$', '$\Delta V_{PVT}$')
ylabel('$\Delta V$')
xlabel('$t$')

figure
hold on
plot(t, mag_p); 
grid on;
ylabel('$\|p\|$')
xlabel('$t$')

figure
hold on
plot(t, log(error)); 
grid on;
ylabel('$e_p$')
xlabel('$t$')

if (parametric)
    figure
    plot(TOF, cost);
    grid on;
    ylabel('cost')
    xlabel('TOF')
end

%% Auxiliary functions 
function [x] = x_update(pInvA, Atb, n, x, z, u) 
   % Impulses update
   x(1:n) = pInvA*(z(1:n)-u(1:n))+Atb;
   dV = reshape(x(1:n), 3, []);
   Norm = sqrt(dot(dV,dV,1));
   dV = -dV./Norm;

   % Primer vector update
   dVm = mean(Norm); 
   sigma = std(Norm);
   eps_V = dVm-2*sigma;

   index = Norm > eps_V;
   for i = 1:size(dV,2)
        if (index(i))
            x(n+1+3*(i-1):n+3*i,1) = dV(:,i);
        else
            x(n+1+3*(i-1):n+3*i,1) = rand(3,1);
        end
   end
end

function [z] = z_update(Phi, b, dVmax, p, n, rho, x, z, u)
    % Impulses update
    start_ind = 1;
    for i = 1:length(p)
        sel = start_ind:p(i);
        
        % L2 minimization
        z(sel) = shrinkage(x(sel) + u(sel), 1/rho);

        % Projection onto the maximum radius balls
        Norm = norm(z(sel));
        if (Norm > dVmax)
           z(sel) = dVmax * z(sel)/norm(z(sel));
        end

        Phi(:,sel) = -Norm * Phi(:,1+3*(i-1):3*i);
        
        start_ind = p(i) + 1;
    end

    % Primer vector updates
    pinvA = pinv(Phi);
    Ab = pinvA*b;
    z(n+1:2*n) = (eye(n)-pinvA*Phi)*(x(n+1:2*n)+u(n+1:2*n))+Ab;
end

% Euclidean norm proximal operator
function z = shrinkage(x, kappa)
    % L2 minimization
    z = max(0, 1 - kappa/norm(x))*x;
end

% Objective
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