%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Primal Solver %% 
% Main solver function %

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
        stm_idx = 1 + m * (i - 1): m * i;
        state_idx = 1 + n * (i - 1) : n * i;

        Phi(:,state_idx) = ( M / STM(:,stm_idx) ) * B(:,state_idx);
    end

    % Compute the initial missvector
    b = xf - M * x0;

    % Pre-factoring of constants
    nx = n * N;
    Id = eye(nx);
    c = zeros(nx,1);

    % Equilibration
    [~, ePhi, ~, D1, D2] = src.RuizEquil( zeros(size(Phi,2),1), Phi, 1E-6, 'L' );
    eb = (D1 .* b.').';

    umax = obj.Actuator.umax;
    umin = obj.Actuator.umin;

    % Normal equations
    invPhi = pinv(ePhi);
    Atb    = invPhi * eb;
    pInvA  = Id - invPhi * ePhi;

    % Create the functions to be solved 
    Obj      = @(x,z)  ( obj.objective(obj.Actuator.p, x) );
    X_update = @(x,z,u)( obj.x_update(pInvA, Atb, x, z, u) );
    Z_update = @(x,z,u)( obj.z_update(n, obj.Actuator.p, obj.Actuator.q, umin, umax, obj.Mission.N, rho, x, z, u) );

    % ADMM consensus constraint definition 
    A = Id;
    B = -A;        

    % Problem
    Problem = src.SolverADMM(Obj, X_update, Z_update, rho, A, B, c);

    if ( ~exist('alpha', 'var') )
        alpha = 1;
    end

    Problem.alpha = alpha;
    Problem.QUIET = false;

    % Optimization
    tic
    [x, ~, Output] = Problem.solver();
    obj.SolveTime = toc;

    % Undo the equilibration 
    x(:,end) = ( D2 .* x(:,end).' ).';

    % Output 
    u = reshape(x(:,end), n, []);        % Control sequence

    obj.Cost = obj.Actuator.p.ComputeVectorNorm( u );
    obj.Cost = sum( obj.Cost );

    obj.Report = Output;                 % Optimization report
    e = b - Phi * x(:,end);              % Final missvector
    obj.e(:,1) = e;                      % Final missvector 
    obj.e(:,2) = obj.e(:,1);             % Final missvector 
    obj.u = u;                           % Final rendezvous impulsive sequence
    obj.t = t;                           % Execution times
end