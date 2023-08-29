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
%          - scalar epsilon, numerical tolerance
%          - scalar rho, the augmented Lagrangian penalty parameter (rho > 0)
%          - scalar alpha, the overfitting parameter (2 > alpha > 0)

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj, epsilon, rho, alpha)
    % Preallocation 
    x0 = obj.Mission.x0;                    % Initial conditions 
    xf = obj.Mission.xf;                    % Final conditions 

    t = obj.Mission.t;                      % Mission clock
    N = length(t);                          % Number of total opportunities
    STM = obj.Mission.Phi;                  % STM of the system
    B = obj.Mission.B;                      % Control input of the system
    m = obj.Mission.m;                      % State vector dimension
    n = obj.Mission.n;                      % Control input dimension

    Phi = zeros(size(B,2), size(STM,1));
    M = STM(:,1+m*(N-1):m*N);

    for i = 1:length(t)
        Phi(1+n*(i-1):n*i,:) = (STM(:,1+m*(i-1):m*i)^(-1)*B(:,1+n*(i-1):n*i)).';
    end

    % Compute the initial missvector
    b = (M\xf) - x0;

    % Save the initial STM
    M = Phi;

    % Optimization of the Lagrange multiplier
    maxIter = 10; 
    iter = 1;
    GoOn = true;

    while (iter < maxIter && GoOn && N >= 1)
        v = [b; zeros(n * N,1)];
    
        % Pre-factoring of constants
        p = repmat(n, 1, N);
        cum_part = cumsum(p);
        pPhi = [Phi kron(eye(N),-eye(n))];
        Theta = [rho * eye(size(pPhi,2)) pPhi.'; pPhi zeros(size(pPhi,1))];
        Theta = pinv(Theta);
    
        % Create the functions to be solved 
        Obj = @(x,z)(obj.objective(v, z));
        X_update = @(x,z,u)(obj.x_update(m, Theta, v, rho, x, z, u));
        Z_update = @(x,z,u)(obj.z_update(cum_part, obj.Thruster.q, obj.Mission.N, Phi, rho, x, z, u));
    
        % ADMM consensus constraint definition 
        A = eye(m + n * N);
        B = -eye(m + n * N);        
        c = zeros(m + n * N,1);
    
        % Problem
        Problem = ADMM_solver(Obj, X_update, Z_update, rho, A, B, c);
    
        if (~exist('alpha', 'var'))
            alpha = 1;
        end
        Problem.alpha = alpha;

        % Solve the problem
        tic
        [x, ~, Output] = Problem.solver();
        obj.SolveTime = toc;

        % Output
        lambda = reshape(x(1:m,end), 1, []).';  % Lagrange multiplier
        p = Phi * lambda;
        p = reshape(p, n, N);

        switch (obj.Thruster.q)
            case 'L2'
                p_norm = sqrt( dot(p,p,1) );
            case 'L1'
                p_norm = sum( abs(p), 1 );
            case 'Linfty'
                p_norm = max( abs(p), [], 1 );
        end

        % Check for convergence
        max_p = sort(p_norm);
        max_p = max_p(end)-1;
        if (max_p > epsilon)
            GoOn = false;
        else
            iter = iter+1;
            index = abs(p_norm-1) <= epsilon;
            N = sum(index);
            index = kron(index, ones(1,n));
            Phi = Phi(logical(index),:);
        end
    end

    % Computation of the control law
    p = M * lambda;
    u = [lambda; p];
    p = reshape(p, n, []);
    N = size(p,2);

    switch (obj.Thruster.q)
        case 'L2'
            p_norm = sqrt( dot(p,p,1) );
        case 'L1'
            p_norm = sum( abs(p), 1 );
        case 'Linfty'
            p_norm = max( abs(p), [], 1 );
    end

    index = abs(p_norm-1) < epsilon;

    if (sum(index) > m)
        dp = abs(p_norm-1) < epsilon;
        [~, pos] = sort(dp, 'descend');
        index = zeros(1,length(p_norm));
        index(pos(1:m)) = ones(1,m);
        index = logical(index);
    end

    id = kron(index, ones(1,n));
    A = M(logical(id),:);
    dV = zeros(n,N);
    if (~isempty(A))
        dV(:,index) = reshape((A.'\b), n, []); 
    end

    % Output
    switch (obj.Thruster.p)
        case 'L2'
            obj.Cost = sqrt( dot(dV,dV,1) );
        case 'L1'
            obj.Cost = sum( abs(dV), 1 );
        case 'Linfty'
            obj.Cost = max( abs(dV), [], 1 );
    end

    obj.Cost = sum(obj.Cost);
   
    obj.Report = Output;                        % Optimization report
    obj.e(:,1) = b - M.' * reshape(dV, [], 1);  % Final missvector   
    obj.u = dV;                                 % Final rendezvous impulsive sequence
    obj.t = t;                                  % Execution times
    e = obj.e;
end