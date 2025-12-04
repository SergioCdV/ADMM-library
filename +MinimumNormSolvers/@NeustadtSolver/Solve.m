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
%          - bool equil_flag, to use or not Ruiz equilibration

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj, epsilon, rho, alpha, init_guess, equil_flag)
    % Sanity checks 
    if ( ~exist('alpha', 'var') )
        alpha = 1;
    end

    if ( ~exist('init_guess', 'var') )
        init_guess = [];
    end

    if ( ~exist('equil_flag', 'var') )
        equil_flag = false;
    end

    % Pre-allocation 
    x0 = obj.Mission.x0;                    % Initial conditions 
    xf = obj.Mission.xf;                    % Final conditions 

    t = obj.Mission.t;                      % Mission clock
    N = length(t);                          % Number of total opportunities
    STM = obj.Mission.Phi;                  % STM of the system
    B = obj.Mission.B;                      % Control input of the system
    m = obj.Mission.m;                      % State vector dimension
    n = obj.Mission.n;                      % Control input dimension

    LambdaIdx = 1 : m;

    Phi = zeros(size(B,2), size(STM,1));    % Pre-allocation of the STM
    Phi0 = STM(:,LambdaIdx);                % Initial STM

    % Compute the STM
    for i = 1:length(t)
        idx = 1 + n * (i - 1) : n * i;
        Phi(idx,:) = ( STM(:,1+m*(i-1):m*i) \ B(:,idx) ).';
    end

    % Compute the initial missvector
    idx = 1 + m * ( N - 1 ) : m * N;
    M = STM(:,idx);
    b = (M \ xf) - (Phi0 \ x0);

    % Constant matrices  << this is for speed in a computation unit with sufficient RAM
    Ones = ones(1,n);       % Vectors of 1
    Id   = eye(n * N);      % Identity matrix of n * N x n * N
    Os   = zeros(n * N);    % Zero matrix of n * N x n * N

    % Initial indices 
    time_mask = logical( Os(1,1:N) );
    time_mask([1 floor(N/2) end]) = true; 
    time_idx = 1 : N; 

    % Number of impulsive opportunities 
    Nopp = sum(time_mask);
    
    % Complete primer vector system 
    KronEye = -Id(1:n*N,1:n*N);      
    pPhi    = [Phi KronEye];   

    % Cost function
    vinit = [-b; -Os(:,1)];

    % Equilibration 
    if ( equil_flag )
        [ev, egPhi, ~, ~, ~] = src.RuizEquil( vinit, pPhi, 1E-6, 'L' );
    else
        ev        = vinit;
        egPhi     = pPhi;
    end

    % Local STM 
    currPhi = PartitionSTM( time_mask, Ones, egPhi(:,LambdaIdx) );

    % Optimization of the Lagrange multiplier
    maxIter = 20;           % Maximum number of iterations
    iter    = 1;            % Current iteration index
    GoOn    = Nopp >= 2;    % Boolean to control convergence

    tic
    while ( iter < maxIter && GoOn && Nopp > 0 )
        % Constants of the iteration 
        dimPrimer = n * Nopp;
        idx = 1 : dimPrimer;
        nx  = m + dimPrimer;
        NxIdx = 1 : nx;                                

        % Linear cost function
        v        = ev(NxIdx);
        b_primer = Os(idx,1);

        if ( iter == 1 )
            % Initial linear system 
            epPhi = PartitionSTM( time_mask, Ones, egPhi );
            idx   = logical([1:m kron(time_mask,Ones)]);
            epPhi = epPhi(:,idx);

            % Initial Cholesky decomposition
            cholPhi = chol(epPhi * epPhi.', "lower");

        else
            % Update initial guess 
            p            = currPhi * lambda;                         
            init_guess.x = [lambda; reshape(p, [], 1)];    
            init_guess.z = init_guess.x;
        end
    
        % Create the functions to be solved 
        Obj      = @(x,z)  ( obj.objective(nx, v, z) );
        X_update = @(x,z,u)( obj.x_update(cholPhi, epPhi, v, -b_primer, rho, x, z, u) );
        Z_update = @(x,z,u)( obj.z_update(m, n, obj.Actuator.q, v(1:m), rho, x, z, u) );
    
        % ADMM consensus constraint definition 
        A = Id(NxIdx,NxIdx);
        B = -A;        
        c = Os(NxIdx,1);
    
        % Problem solve
        Solv = src.SolverADMM(Obj, X_update, Z_update, rho, A, B, c, init_guess);

        Solv.alpha = alpha;
        Solv.QUIET = false;

        % Solve the problem
        [x, ~, Output] = Solv.solver();

        % Output
        lambda = reshape(x(LambdaIdx,end), 1, []).';    % Lagrange multiplier
        p = egPhi(:,LambdaIdx) * lambda;                % Primer vector
        p = reshape(p, n, N);                           % Primer vector
        
        % Check for convergence
        p_norm = obj.Actuator.q.ComputeVectorNorm( p ); % Switching surface
        [max_p, pos] = sort(p_norm);

        if ( max_p(end) <= 1 + epsilon )
            % Convergence
            GoOn = false;
        else
            % Include the new maximum 
            old_mask = time_mask;
            time_mask( pos(end) ) = true;

            % Do not include the non-plausible actuation epochs
            index = p_norm < 1 - epsilon;                
            time_mask( index ) = 0;

            % Complete matrix
            currPhi = PartitionSTM( time_mask, Ones, egPhi(:,LambdaIdx) );

            % Downdate the STM
            rem_pos = old_mask & ~time_mask;
            rem_pos = time_idx( rem_pos );
            old_pos = time_idx( old_mask );
            Nrm = length( rem_pos );

            if ( Nrm > 0 )
                rem_pos = find( ismember( old_pos, rem_pos ), Nrm );
                rem_pos = (rem_pos-1) * n + (1:n).'; 
                rem_pos = rem_pos(:).';

                % Update the Cholesky factor
                cholPhi = src.RemdateChol( cholPhi, rem_pos );

                epPhi(rem_pos,:) = [];              % Delete rows 
                rem_pos          = m + rem_pos;     % Column indices
                epPhi(:,rem_pos) = [];              % Delete columns

                Nrm = n * Nrm;                      % Number of removed variables
            end

            % Update STM 
            idx     = time_mask & ~old_mask;
            Npls    = sum(idx);
            newPhi  = PartitionSTM( idx, Ones, egPhi(:,LambdaIdx) );

            idx     = 1 : n * Npls;
            newPhi  = [newPhi Os(idx,1:nx-m-Nrm) -Id(idx,idx)];
            epPhi   = [epPhi Os(1:nx-m-Nrm,idx)];

            % Update of the Cholesky decomposition of the STM inverse
            [cholPhi, epPhi] = src.AggdateChol( cholPhi, epPhi, newPhi );
            
            % Number of impulsive opportunities 
            Nopp = sum(time_mask);

            % Update the iteration counter
            iter = iter + 1;
        end
    end

    % Computation of the control law
    if ( 1) %~GoOn )
        % Final output 
        u = [lambda; egPhi(:,LambdaIdx) * lambda];                 % Adjoint vector at final epoch and primer vector

        % Input reconstruction 
        [t_pruned, dv] = obj.ImpulseReconstruction(t, b, Phi, p_norm, 1, epsilon);

        % Complete action sequence
        dV = zeros( n, length(t) );               
        for i = 1:length(t_pruned)
            dV(:, t_pruned(i) == t) = dv(:,i);
        end
    
        % Output
        e = b - Phi.' * reshape(dV, [], 1);         % Regulation error
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

        u  = zeros(m + n * N,1); 
        e  = b;
    end
    obj.SolveTime = toc;
end

%% Auxiliary functions 
function [A] = PartitionSTM(mask, ones, Phi)
    index = kron(mask, ones);                    % Actuation epochs
    A     = Phi(logical(index),:);               % STM corresponding to the new actuation grid
end