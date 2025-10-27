%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Carter Solver %% 
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
        cntrl_idx = 1 + n * (i - 1) :n * i;
        state_idx = 1 + m * (i - 1) : m * i;

        Phi(:,cntrl_idx) = ( STM(:,state_idx) \ M ) * B(:,cntrl_idx);
    end

    % Compute the initial missvector
    b = xf - M * x0;

    % Pre-factoring of constants
    Atb = pinv(Phi) * b;
    pInvA = eye(size(Phi,2)) - pinv(Phi) * Phi;

    % Create the functions to be solved 
    Obj = @(x,z)( obj.objective(obj.Actuator.p, x) );
    X_update = @(x,z,u)( obj.x_update(n, obj.Actuator.q, pInvA, Atb, x, z, u) );
    Z_update = @(x,z,u)( obj.z_update(n, obj.Actuator.p, obj.Actuator.q, obj.Actuator.umin, obj.Actuator.umax, obj.Mission.N, Phi, b, rho, x, z, u) );

    % ADMM consensus constraint definition 
    nx = 2 * n * N;
    A = eye(nx);
    B = -A;        
    c = zeros(nx,1);

    % Problem
    Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);

    if ( ~exist('alpha', 'var') )
        alpha = 1;
    end

    Problem.alpha = alpha;
    Problem.QUIET = false;

    % Optimization
    tic
    [x, z, Output] = Problem.solver();
    obj.SolveTime = toc;

    % Output 
    nv = n * N;
    dV = reshape( x(1:nv,end), n, [] );       % Control sequence
    p  = reshape( x(nv+1:end,end), n, [] );   % Primer vector
    u  = [dV; p]; 

    obj.Cost = obj.Actuator.p.ComputeVectorNorm( u );
    obj.Cost = sum(obj.Cost);

    obj.e(:,1) = b - Phi * x(1:nv,end);       % Final missvector  
    obj.e(:,2) = b - Phi * z(1:nv,end);       % Final missvector 
    e = obj.e;

    obj.Report = Output;                      % Optimization report
    obj.u = dV;                               % Final rendezvous impulsive sequence
    obj.t = t;                                % Execution times
end