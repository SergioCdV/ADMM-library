%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Neustadt Solver %% 
% Main solver function %

% Inputs:  - object obj, the Linear Rendezvous Problem object
%          - scalar epsilon, numerical tolerance
%          - scalar rho, the augmented Lagrangian penalty parameter (rho > 0)
%          - scalar alpha, the overfitting parameter (2 > alpha > 0)
%          - vector init_guess, an initial guess for the optimization

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj, epsilon, rho, alpha, init_guess)
    % Sanity checks 
    if ( ~exist('alpha', 'var') )
        alpha = 1;
    end

    % Pre-allocation 
    x0 = obj.Mission.x0;                    % Initial conditions 
    xf = obj.Mission.xf;                    % Final conditions 

    t = obj.Mission.t;                      % Mission clock
    t_pruned = t;                           % Mission clock
    N = length(t);                          % Number of total opportunities
    STM = obj.Mission.Phi;                  % STM of the system
    B = obj.Mission.B;                      % Control input of the system
    m = obj.Mission.m;                      % State vector dimension
    n = obj.Mission.n;                      % Control input dimension

    Phi = zeros(size(B,2), size(STM,1));    % Pre-allocation of the STM
    Phi0 = STM(:,1:m);                      % Initial STM

    % Compute the STM
    for i = 1:length(t)
        idx = 1+n*(i-1):n*i;
        Phi(idx,:) = ( STM(:,1+m*(i-1):m*i) \ B(:,idx) ).';
    end

    % Compute the initial missvector
    idx = 1 + m * ( N - 1 ) : m * N;
    M = STM(:,idx);
    b = (M \ xf) - (Phi0 \ x0);

    % Constant matrices  << this is for speed in a computation unit with sufficient RAM
    M    = Phi;             % Initial STM
    Ones = ones(1,n);       % Vectors of 1
    Id   = eye(n);          % Identity matrix of n x n
    Os   = zeros(n * N);    % Zero matrix of n * N x n * N

    % Optimization of the Lagrange multiplier
    maxIter = 20;           % Maximum number of iterations
    iter    = 1;            % Current iteration index
    GoOn    = N >= 2;       % Boolean to control convergence

    % Initial indices 
    time_mask = logical( Os(1,1:N) );
    time_mask(1) = true; 
    time_mask(end) = true; 
    time_mask( floor(N/2) ) = true;

    % Cost function
    vinit = [-b; -Os(:,1)];

    while ( iter < maxIter && GoOn )
        % Number of impulsive opportunities 
        Nopp = sum(time_mask);

        % Local STM 
        index    = kron(time_mask, Ones);               % Actuation epochs
        curr_Phi = Phi(logical(index),:);               % STM corresponding to the new actuation grid

        % Primer vector linear system
        KronEye = kron( eye(Nopp), -Id );                            
        pPhi = [curr_Phi KronEye];                                    
        idx = 1 : n * Nopp;
        Theta = [rho * eye(size(pPhi,2)) pPhi.'; pPhi Os(idx,idx)];
        Theta = pinv(Theta);

        % Linear cost function at each iteration grid
        nx = m + n * Nopp;
        v = vinit(1:nx);
        linear_cost = [v; -zeros(size(Theta,1)-nx,1)];
    
        % Create the functions to be solved 
        Obj = @(x,z)( obj.objective(v, z) );
        X_update = @(x,z,u)( obj.x_update( Theta, linear_cost, rho, x, z, u) );
        Z_update = @(x,z,u)( obj.z_update( n, obj.Actuator.q, Phi, -b, rho, x, z, u) );
    
        % ADMM consensus constraint definition 
        A = eye(nx);
        B = -A;        
        c = zeros(nx,1);
    
        % Problem
        if ( iter == 1 && exist( 'init_guess', 'var' ) )
            Solv = src.SolverADMM(Obj, X_update, Z_update, rho, A, B, c, init_guess);
        else
            Solv = src.SolverADMM(Obj, X_update, Z_update, rho, A, B, c);
        end

        Solv.alpha = alpha;
        Solv.QUIET = false;

        % Solve the problem
        tic
        [x, ~, Output] = Solv.solver();
        obj.SolveTime = toc;

        % Output
        lambda = reshape(x(1:m,end), 1, []).';          % Lagrange multiplier
        p = Phi * lambda;                               % Primer vector
        p = reshape(p, n, N);                           % Primer vector
        
        % Check for convergence
        p_norm = obj.Actuator.q.ComputeVectorNorm( p ); % Switching surface
        [max_p, pos] = sort(p_norm);

        if ( max_p(end) <= 1 + epsilon )
            % Convergence
            GoOn = false;
        else
            % Include the new maximum 
            time_mask( pos(end) ) = true;

            % Do not include the non-plausible actuation epochs
            index = p_norm < 1 - epsilon;                
            time_mask( index ) = zeros(1, sum(index));

            % Update the iteration counter
            iter = iter + 1;
        end
    end

    % Computation of the control law
    if ( Output.Result )
        % Final output 
        u = [lambda; M * lambda];       % Adjoint vector at final epoch and primer vector

        % Input reconstruction 
        [t_pruned, dv] = obj.ImpulseReconstruction(t, b, Phi, p_norm, epsilon);

        % Complete action sequence
        dV = zeros( n, length(t) );               
        for i = 1:length(t_pruned)
            dV(:, t_pruned(i) == t) = dv(:,i);
        end
    
        % Output
        e = b - M.' * reshape(dV, [], 1);           % Regulation error
        obj.e(:,1) = e;                             % Final missvector   
        obj.Cost = dot(b, lambda);                  % Final minimum-norm cost
        obj.Report = Output;                        % Optimization report
        obj.u = dV;                                 % Final impulsive sequence
        obj.t = t;                                  % Execution times
        
    else
        % TODO: see what is needed here
        % Check the q-norm of the primer vector on the z-update
%         p = reshape( z(m+1:end,end), n, [] );
%         p_norm = obj.Actuator.q.ComputeVectorNorm( p );

        u  = []; 
        e  = [];
    end
end