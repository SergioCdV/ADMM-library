%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Hybrid Solver %% 
% Main solver function %

% Inputs:  - object obj, the Linear Rendezvous Problem object
%          - scalar rho, the augmented Lagrangian penalty parameter (rho > 0)
%          - scalar alpha, the overfitting parameter (2 > alpha > 0)

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj, rho, alpha, equil_flag)
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
        cntrl_index = 1 + n * (i - 1) : n * i;
        state_idx = 1 + m * (i - 1) : m * i;

        Phi(:,cntrl_index) = ( M / STM(:,state_idx) ) * B(:,cntrl_index);
    end

    % Compute the initial missvector
    b = xf - M * x0;

    % Pre-factoring of constants
    nx = n * N;
    Id = eye(nx);
    c = zeros(nx,1);

    A = [(1 + rho) * eye(size(Phi,2)) Phi.'; Phi zeros(size(Phi,1))];
    A = pinv(A);

    % Equilibration
    if ( ~exist('equil_flag', 'var') )
        equil_flag = true;
    end

    if ( equil_flag )
        [~, eA, ~, D1, ~] = src.RuizEquil( zeros(size(Phi,2),1), A, 1E-6, 'L' );
        eb = (D1 .* b.').';
    else
        eA = A; 
        eb = b;
    end

    umax = obj.Actuator.umax;
    umin = obj.Actuator.umin;

    % Create the functions to be solved 
    Obj = @(x,z)( obj.objective(x, z) );
    X_update = @(x,z,u)( obj.x_update(eA, eb, rho, x, z, u) );
    Z_update = @(x,z,u)( obj.z_update(n, umin, umax, obj.Mission.N, rho, x, z, u) );

    % ADMM consensus constraint definition 
    A = Id;
    B = -A;        

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
    dV2 = reshape(x(:,end), n, []);      % Control sequence
    dV  = reshape(z(:,end), n, []);      % Control sequence
    u   = [dV2; dV]; 

    obj.Cost = sum( abs(dV), 1 );
    obj.Cost = sum( obj.Cost );

    obj.e(:,1) = b - Phi * x(:,end);     % Final missvector  
    obj.e(:,2) = b - Phi * z(:,end);     % Final missvector 
    e = obj.e;

    obj.Report = Output;                 % Optimization report
    obj.u = dV;                          % Final rendezvous impulsive sequence
    obj.t = t;                           % Execution times
end