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

    umax = obj.Actuator.umax;                   % Maximum control authority

    if ( umax == Inf )
        % Call the standard Neustadt solver 
        inf_solver = MinimumNormSolvers.NeustadtSolver(obj.Mission, obj.Actuator);
        [t, u, e, obj] = inf_solver.Solve(epsilon, rho, alpha, init_guess);
        
    else    
        % Pre-allocation 
        x0 = obj.Mission.x0;                    % Initial conditions 
        xf = obj.Mission.xf;                    % Final conditions 
    
        t = obj.Mission.t;                      % Mission clock
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

        % Initial time window and the corresponding indices
        w = obj.window_ratio * (t(end) - t(1));                 % Time duration of each window
        Nw = floor( 1 / obj.window_ratio ) - 1;                 % Window index
        edges = [t(1) + w * (0:Nw) t(end)];                     % Window edges
        tw_idx = discretize(t, edges);                          % Map between t and each window
        
        [unique_idx, first_idx] = unique(tw_idx, 'first');      % Index of the start of each window    
        [~, last_idx]           = unique(tw_idx, 'last');       % Index of the end of each window
        
        edge_idx = [first_idx, last_idx];                       % Pair the indices of each window
        edge_idx = reshape( edge_idx.', 1, [] );                % Reshaping of start-end indices

        Tk = numel(unique_idx);                                 % Number of active windows
        Sigma = zeros(1,Tk);                                    % Original vector of Lagrange multipliers associated to the bound constraint
        Slack = zeros(1,Tk);                                    % Original vector of slack Lagrange multipliers associated to the bound constraint
        
        % Constant matrices  << this is for speed in a computation unit with sufficient RAM
        Nx = m + 2 * Tk + n * N;                                % Original number of variables
        Ones = ones(1,Nx);                                      % Vectors of 1
        Id   = eye(Nx);                                         % Identity matrix of n x n
        Os   = zeros(Nx);                                       % Zero matrix of nN x nN

        % Initial time mask
        time_mask = logical( Os(1,1:N) );                       % Pre-allocation      
        time_mask(edge_idx) = true;                             % Initial time mask
        Nopp = sum(time_mask);                                  % Number of impulsive opportunities 

        % Local STM                 
        index    = kron(time_mask, Ones(1,1:n));                % Actuation epochs
        curr_Phi = Phi(logical(index),:);                       % STM corresponding to the new actuation grid  

        % Complete cost function
        vinit = [-b; +umax * ones(Tk,1); zeros(Tk,1); +Os(:,1)]; 
        
        % Index mapping of variables
        lambda_pos = 1 : m;                                     % Position of the Lagrange multiplier
        sigma_pos  = m + 1 : m + Tk;                            % Position of the bound Lagrange multipliers
        primer_pos = Nx + 1 - n * Nopp : Nx;                    % Position of the primer vector
        
        sigma_map  = tw_idx(time_mask);                         % Mapping between primer vector and Lagrange multipliers
        [~, sigma_unique] = unique( sigma_map, 'first' );       % Mapping between primer vector and unique Lagrange multipliers

        % Optimization of the Lagrange multiplier
        maxIter = 10;                                           % Maximum number of iterations
        iter    = 1;                                            % Current iteration index
        GoOn    = N >= 2;                                       % Boolean to control convergence

        while ( iter < maxIter && GoOn && Nopp > 0 )
            % Primer vector linear system
            KronEye = kron( eye(Nopp), -Id(1:n,1:n) );  
            idx = 1 : n * Nopp;
            KronZero = Os(idx,1: 2 * Tk);
            primer_system = [curr_Phi KronZero KronEye];

            % Augmented slack variables system 
            slack_system = [Os(1:Tk,1:m) -Id(1:Tk,1:Tk) Id(1:Tk,1:Tk) Os(1:Tk,idx)];
            
            % Pre-allocation of the pseudoinverse of the linear inverse
            pPhi = [primer_system; slack_system];
            Theta = [rho * eye(size(pPhi,2)) pPhi.'; pPhi zeros(size(pPhi,1))];
            Theta = pinv(Theta);
    
            % Linear cost function at each iteration grid
            nx = m + 2 * Tk + n * Nopp;                                            % Number of decision variables
            v = vinit( [lambda_pos sigma_pos Tk + sigma_pos primer_pos] );         % Current cost function
            linear_cost = [v; -zeros(size(Theta,1)-nx-Tk,1); -Ones(1,1:Tk).'];     % KKT cost function
    
            % Create the functions to be solved 
            Obj = @(x,z)( MinimumNormSolvers.NeustadtSolver.objective(v, z) );
            X_update = @(x,z,u)( obj.x_update( Theta, linear_cost, sigma_pos, rho, x, z, u ) );
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
            lambda = x(lambda_pos,end);          % Lagrange multiplier
            sigma  = x(sigma_pos ,end);          % Lagrange multiplier associated to the control bound
            slackT = x(Tk + sigma_pos,end);      % Slack variables associated to the Lagrange multipliers
            p      = Phi * lambda;               % Primer vector
            p      = reshape(p, n, N);           % Primer vector

            % Update the complete set of Lagrange multipliers
            Sigma(unique_idx) = sigma;  
            Slack(unique_idx) = slackT;
            
            % Check for convergence
            p_norm = obj.Actuator.q.ComputeVectorNorm( p );        % Switching surface
          
            if ( all( p_norm <= (1 + Sigma(tw_idx)) + epsilon ) && Output.Result )
                % Convergence
                GoOn = false;
            else
                % Time-window analysis
                for i = 1:length(Sigma)
                    % Include the new maximum per time window
                    pos = 1 + 2 * (i - 1) : 2 * i;
                    range = edge_idx(pos);
                    tw_norm_idx = range(1) : range(2);
                    tw_norm = p_norm( tw_norm_idx );
                    [~, pos] = sort( tw_norm );
                    time_mask( tw_norm_idx( pos(end) ) ) = true;

                    % Do not include the non-plausible actuation epochs
%                     discard_idx = tw_norm < 1 - epsilon;                
%                     time_mask( tw_norm_idx(discard_idx) ) = 0;
                end
    
                % Update the time window
                [unique_idx, sigma_unique, sigma_map] = unique( tw_idx( time_mask ), 'first' );
                Tk = numel(unique_idx);

                % Update the number of variables
                Nopp = sum(time_mask);                              % Number of impulsive opportunities 
                sigma_pos  = m + 1 : m + Tk;                        % Position of the bound Lagrange multipliers
                primer_pos = Nx + 1 - n * Nopp : Nx;                % Position of the primer vector
            
                % Local STM 
                index    = kron(time_mask, Ones(1,1:n));            % Actuation epochs
                curr_Phi = Phi(logical(index),:);                   % STM corresponding to the new actuation grid

                % Update initial guess
                p = curr_Phi * lambda;                              % Initial guess for the primer vector
                sigma = Sigma(unique_idx);                          % Initial guess for the Lagrange multipliers
                slackT = Slack(unique_idx);                         % Initial guess for the slack variables

                init_guess.x = [lambda; sigma.'; slackT.'; reshape(p, [], 1)];
                init_guess.z = init_guess.x;

                init_guess = [];

                % Update the iteration counter
                iter = iter + 1;
            end
        end
        
        % Computation of the control law
        if ( 1 )%~GoOn )
            % Final output 
            u = [lambda; Phi * lambda];       % Adjoint vector at final epoch and primer vector
    
            % Input reconstruction 
            [t_pruned, dv] = MinimumNormSolvers.NeustadtSolver.ImpulseReconstruction(t, b, Phi, p_norm, 1 + Sigma(tw_idx), epsilon);
    
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
    end
end