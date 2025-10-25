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
    p = repmat(n, 1, N);    % Control indices across time
    cum_part = cumsum(p);   % Control indices across time

    % Optimization of the Lagrange multiplier
    maxIter = 10;           % Maximum number of iterations
    iter    = 1;            % Current iteration index
    GoOn    = true;         % Boolean to control convergence

    % Cost function
    vinit = [-b; -zeros(n * N,1)];

    while (iter < maxIter && GoOn && N > 1)
        % Primer vector linear system
        KronEye = kron(eye(N), -Id);                            
        pPhi = [Phi KronEye];                                                   
        Theta = [rho * eye(size(pPhi,2)) pPhi.'; pPhi Os(1:n*N, 1:n*N)];
        Theta = pinv(Theta);

        % Linear cost function at each iteration grid
        nx = m + n * N;
        v = vinit(1:nx);
        linear_cost = [v; -zeros(size(Theta,1)-nx,1)];
    
        % Create the functions to be solved 
        Obj = @(x,z)( obj.objective(v, z) );
        X_update = @(x,z,u)( obj.x_update( Theta, linear_cost, rho, x, z, u) );
        Z_update = @(x,z,u)( obj.z_update( cum_part(1:N), obj.Actuator.q, Phi, -b, rho, x, z, u) );
    
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
        p_norm = obj.Actuator.q.ComputeVectorNorm( p ); % Switching surface
        
        % Check for convergence
        max_p = sort(p_norm);

        if ( max_p(end) <= 1 + epsilon(1) )
            % Convergence
            GoOn = false;
        else
            % Reduce the number of points to analyze
            index = p_norm >= 1 - epsilon(2);           % Plausible actuation epochs
            N = sum(index);                             % Density of the grid
            t_pruned = t_pruned( logical(index)  );     % Actuation epochs
            index = kron(index, Ones);                  % Actuation epochs
            Phi = Phi(logical(index),:);                % STM corresponding to the new actuation grid

            % Update the iteration counter
            iter = iter + 1;
        end
    end
    
    % Final output 
    u = [lambda; M * lambda];       % Adjoint vector at final epoch and primer vector

    % Computation of the control law
    if ( Output.Result )
        % Both the z and x solutions are equivalent
        N       = size(p,2);
        imp_opp = abs(p_norm - 1) <= epsilon(1);

        if ( sum(imp_opp) > 1E3 )
     % TODO: analysis based on the derivative of the primer vector
    %         dp = p_norm - 1;
    %         d_dp  = gradient(dp); 
    %         dd_dp = gradient(d_dp);
    %         dp = abs(dp);
    %         [~, pos] = sort(dp);
    %         index = zeros(1,length(p_norm));
    % 
    %         k = 0;
    %         for i = 1:length(pos)
    %             if (t(pos(i)) == t(1) && sign(d_dp(1)) < 0)
    %                 index(pos(i)) = 1;
    %                 k = k+1;
    %             elseif (t(pos(i)) == t(end) && sign(d_dp(end)) > 0)
    %                 index(pos(i)) = 1;
    %                 k = k+1;
    %             elseif (sign(dd_dp(pos(i))) < 0)
    %                 index(pos(i)) = 1;
    %                 k = k+1;
    %             end
    % 
    %             if (k == m)
    %                 break;
    %             end
    %         end
    %         
    %         index = logical(index);
    %         t_pruned = t_pruned(index);
    %         index = kron(index, ones(1,n));

        elseif ( ~isempty(Phi) )
            % Check actuation epochs
            t_pruned = t_pruned( logical(imp_opp) );
            index    = kron( imp_opp, Ones );
            index    = logical( index );
        end
    
        % Action sequence
        if ( ~isempty(Phi) )
            % Compute the maneuver sequence via solving the corresponding linear system
            dv = Phi(index,:).' \ b; 
            dv = reshape( dv, n, [] );
            
            % Complete action sequence
            dV = zeros(n, length(t));               
            for i = 1:length(t_pruned)
                dV(:, t_pruned(i) == t) = dv(:,i);
            end
            
        else
            dV = zeros(n, N);
        end
    
        % Output
        e = b - M.' * reshape(dV, [], 1);           % Regulation error
        obj.Cost = dot(b, lambda);                  % Final minimum-norm cost
        obj.Report = Output;                        % Optimization report
        obj.e(:,1) = e;                             % Final missvector   
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