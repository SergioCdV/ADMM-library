%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Dual control-bounded Solver %% 
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

    if ( ~exist('init_guess', 'var') )
        init_guess = [];
    end

    umax = obj.Actuator.umax;               % Maximum control authority

    if ( umax == Inf )
        % Call the standard Neustadt solver 
        inf_solver = MinimumNormSolvers.NeustadtSolver(obj.Mission, obj.Actuator);
        [t, u, e, obj] = inf_solver.Solve(epsilon, rho, alpha, init_guess);
        
    else    
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
        Os   = zeros(n * N);    % Zero matrix of m x m

        % Initial indices and time windows
        w = obj.window_ratio * (t(end) - t(1));
        Nw = floor( 1 / obj.window_ratio ) - 1;
        edges = [t(1) + w * (0:Nw) t(end)];
        tw_idx = discretize(t, edges);
        [unique_idx, first_idx] = unique(tw_idx, 'first');
        [~, last_idx]           = unique(tw_idx, 'last');
        indices = [first_idx, last_idx];
        indices = reshape( indices.', 1, [] );
        Tk = max(unique_idx);

        time_mask = logical( Os(1,1:N) );
        time_mask(indices) = true * ones(1,length(indices)); 

        % Number of impulsive opportunities 
        Nopp = sum(time_mask);
            
        % Local STM 
        index    = kron(time_mask, Ones);                       % Actuation epochs
        curr_Phi = Phi(logical(index),:);                       % STM corresponding to the new actuation grid  

        % Cost function
        vinit = [-b; umax * 0 * ones(N,1); -Os(:,1)];           % Complete cost function
        
        % Index mapping of variables
        Nx = m + Tk + n * N;                                    % Original number of variables
        lambda_pos = 1 : m;                                     % Position of the Lagrange multiplier
        sigma_pos  = m + 1 : m + Tk;                            % Position of the bound Lagrange multipliers
        primer_pos = Nx + 1 - n * Nopp : Nx;                    % Position of the primer vector
        sigma_map = tw_idx(time_mask);                          % Mapping between primer vector and Lagrange multipliers
        [~, sigma_unique] = unique( sigma_map, 'first' );       % Mapping between primer vector and unique Lagrange multipliers
        Sigma = zeros(Tk, 1);                                   % Original vector of Lagrange multipliers associated to the bound constraint

        % Optimization of the Lagrange multiplier
        maxIter = 10;           % Maximum number of iterations
        iter    = 1;            % Current iteration index
        GoOn    = N >= 2;       % Boolean to control convergence

        while ( iter < maxIter && GoOn )
            % Primer vector linear system
            KronEye = kron( eye(Nopp), -Id );  
            idx = 1 : n * Nopp;
            KronZero = Os(idx,1:Tk);
            pPhi = [curr_Phi KronZero KronEye];                                                   
            Theta = [rho * eye(size(pPhi,2)) pPhi.'; pPhi Os(idx,idx)];
            Theta = pinv(Theta);
    
            % Linear cost function at each iteration grid
            nx = m + Tk + n * Nopp;
            v = vinit( [lambda_pos sigma_pos primer_pos] );         % Initial cost function
            linear_cost = [v; -zeros(size(Theta,1)-nx,1)]; 
    
            % Create the functions to be solved 
            Obj = @(x,z)( MinimumNormSolvers.NeustadtSolver.objective(v, z) );
            X_update = @(x,z,u)( MinimumNormSolvers.NeustadtSolver.x_update( Theta, linear_cost, rho, x, z, u ) );
            Z_update = @(x,z,u)( obj.z_update( n, sigma_pos, sigma_map, sigma_unique, obj.Actuator.q, Phi, -b, rho, x, z, u ) );
        
            % ADMM consensus constraint definition 
            A = eye(nx);
            B = -A;        
            c = zeros(nx,1);
        
            % Problem solve
            Solv = src.SolverADMM(Obj, X_update, Z_update, rho, A, B, c, init_guess);
    
            Solv.alpha = alpha;
            Solv.QUIET = false;
    
            % Solve the problem
            tic
            [x, ~, Output] = Solv.solver();
            obj.SolveTime = toc;
    
            % Output
            lambda = reshape(x(lambda_pos,end), 1, []).';          % Lagrange multiplier
            
            sigma = reshape(x(sigma_pos,end), 1, []).';      % Lagrange multiplier associated to the control bound
%             Sigma(sigma_) = sigma;                       % Update the complete set of Lagrange multipliers
            
            p = Phi * lambda;                               % Primer vector
            p = reshape(p, n, N);                           % Primer vector
            
            % Check for convergence
            p_norm = obj.Actuator.q.ComputeVectorNorm( p ); % Switching surface
          
            if ( any( p_norm > (1 + sigma(sigma_map)) + epsilon ) )
                % Time-window analysis
                for i = 1:Tk
                    % Include the new maximum per time window
                    pos = 1 + 2 * (i-1):2 * i;
                    range = indices(pos);
                    tw_norm = p_norm( range(1):range(2) );
                    [~, pos] = sort( tw_norm );
                    time_mask( pos(end) ) = true;

                    % Do not include the non-plausible actuation epochs
                    index = tw_norm < (1 + sigma(i)) - epsilon;                
                    time_mask( index ) = zeros(1, sum(index));
                end
    
                % Update the time window
                sigma_map = tw_idx(time_mask);
                [unique_idx, sigma_unique] = unique(sigma_map, 'first');
                Tk = max(unique_idx);

                Nopp = sum(time_mask);                              % Number of impulsive opportunities 
                sigma_pos  = m + 1 : m + Tk;                        % Position of the bound Lagrange multipliers
                primer_pos = Nx + 1 - n * Nopp : Nx;                % Position of the primer vector
            
                % Local STM 
                index    = kron(time_mask, Ones);                   % Actuation epochs
                curr_Phi = Phi(logical(index),:);                   % STM corresponding to the new actuation grid

                % Update initial guess
                p = curr_Phi * lambda;                              % Initial guess for the primer vector
                sigma = Sigma(sigma_pos);                           % Initial guess for the Lagrange multipliers

                init_guess.x = [lambda; sigma; reshape(p, [], 1)];
                init_guess.x = init_guess.z;

                % Update the iteration counter
                iter = iter + 1;
            else
                % Convergence
                GoOn = false;
            end
        end
        
        % Computation of the control law
        if ( Output.Result )
            % Final output 
            u = [lambda; M * lambda];       % Adjoint vector at final epoch and primer vector
    
            % Input reconstruction 
            [t_pruned, dv] = MinimumNormSolvers.NeustadtSolver.ImpulseReconstruction(t, b, Phi, p_norm, epsilon);
    
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
end