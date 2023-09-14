%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Primal Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

% Inputs:  - object obj, the Linear Rendezvous Problem object
%          - scalar rho, the augmented Lagrangian penalty parameter (rho > 0)
%          - scalar alpha, the overfitting parameter (2 > alpha > 0)

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj, rho, alpha)
    % Preallocation 
    x0 = obj.Mission.x0;                    % Initial conditions 
    xf = obj.Mission.xf;                    % Final conditions 

    t = obj.Mission.t;                      % Mission clock
    N = length(t);                          % Number of total opportunities
    STM = obj.Mission.Phi;                  % STM of the system
    B = obj.Mission.B;                      % Control input of the system
    m = obj.Mission.m;                      % State vector dimension
    n = obj.Mission.n;                      % Control input dimension

    Phi = zeros(size(STM,1), size(B,2));
    M = STM(:,1+m*(N-1):m*N);

    for i = 1:length(t)
        Phi(:,1+n*(i-1):n*i) = (STM(:,1+m*(i-1):m*i)\M) * B(:,1+n*(i-1):n*i);
    end

    % Compute the initial missvector
    b = xf - M * x0;

    % Pre-factoring of constants
    Atb = pinv(Phi)*b;
    pInvA = (eye(size(Phi,2))-pinv(Phi)*Phi);
    p = repmat(n, 1, N);
    cum_part = cumsum(p);

    % Create the functions to be solved 
    Obj = @(x,z)(obj.objective(obj.Thruster.p, cum_part, x));
    X_update = @(x,z,u)(obj.x_update(n, obj.Thruster.q, pInvA, Atb, x, z, u));
    Z_update = @(x,z,u)(obj.z_update(cum_part, obj.Thruster.p, obj.Thruster.q, obj.Thruster.umin, obj.Thruster.umax, obj.Mission.N, Phi, b, rho, x, z, u));

    % ADMM consensus constraint definition 
    A = eye(2 * n * N);
    B = -eye(2 * n * N);        
    c = zeros(2 * n * N,1);

    % Problem
    Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);

    if (~exist('alpha', 'var'))
        alpha = 1;
    end
    Problem.alpha = alpha;
    Problem.QUIET = false;

    % Optimization
    tic
    [x, z, Output] = Problem.solver();
    obj.SolveTime = toc;

    % Output 
    dV = reshape(x(1:n*N,end), n, []);      % Control sequence
    p = reshape(x(n*N+1:end,end), n, []);   % Primer vector
    u = [dV; p]; 

    switch (obj.Thruster.p)
        case 'L2'
            obj.Cost = sqrt( dot(dV,dV,1) );
        case 'L1'
            obj.Cost = sum( abs(dV), 1 );
        case 'Linfty'
            obj.Cost = max( abs(dV), [], 1 );
    end

    obj.Cost = sum(obj.Cost);

    obj.Report = Output;                     % Optimization report
    obj.e(:,1) = b - Phi * x(1:n*N,end);     % Final missvector  
    obj.e(:,2) = b - Phi * z(1:n*N,end);     % Final missvector 
    obj.u = dV;                              % Final rendezvous impulsive sequence
    obj.t = t;                               % Execution times
    e = obj.e;
end