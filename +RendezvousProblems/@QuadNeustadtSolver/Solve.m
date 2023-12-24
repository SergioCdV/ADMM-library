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
    Phi0 = STM(:,1:m);
    M = STM(:,1+m*(N-1):m*N);

    for i = 1:length(t)
        Phi(1+n*(i-1):n*i,:) = (STM(:,1+m*(i-1):m*i)\B(:,1+n*(i-1):n*i)).';
    end

    % Independent grid 
    t = obj.Mission.t;
    t_pruned = t;

    % Compute the initial missvector
    b = (M\xf) - (Phi0\x0);

    % Save the initial STM
    M = Phi;

    % Optimization of the Lagrange multiplier
    maxIter = 10; 
    iter = 1;
    GoOn = true;

    while (iter < maxIter && GoOn && N > 1)
        v = -b;
    
        % Pre-factoring of constants
        p = repmat(n, 1, N);
        cum_part = cumsum(p);
        pPhi = [-b.'; Phi];
    
        % Create the functions to be solved 
        CC = @(z)(obj.ConeCheck(cum_part, obj.Thruster.q, -b, z));
        CP = @(z)(obj.ConeProj(cum_part, obj.Thruster.q, b, z));
        
        % Problem
        Problem = Solvers.AQP_solver(CC, CP, zeros(size(Phi,2), size(Phi,2)), -b, pPhi);
    
        if (~exist('alpha', 'var'))
            alpha = 1;
        end
        Problem.alpha = alpha;
        Problem.QUIET = false;

        % Solve the problem
        tic
        [x, z, Output] = Problem.solver();
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
        if (max_p(end) <= 1 + epsilon(1))
            GoOn = false;
        else
            iter = iter+1;
            index = p_norm > 1 - epsilon(2);
            N = sum(index);
            t_pruned = t_pruned(logical(index));
            index = kron(index, ones(1,n));
            Phi = Phi(logical(index),:);
        end
    end
    
    % Final output 
    u = [lambda; M * lambda];

    % Computation of the control law
    p = reshape(z(m+1:end,end), n, []);
    N = size(p,2);

    switch (obj.Thruster.q)
        case 'L2'
            p_norm = sqrt( dot(p,p,1) );
        case 'L1'
            p_norm = sum( abs(p), 1 );
        case 'Linfty'
            p_norm = max( abs(p), [], 1 );
    end

    imp_opp = abs(p_norm-1) < epsilon(1);

    if (sum(imp_opp) > 1e3)
        dp = p_norm-1;
        d_dp = gradient(dp); 
        dd_dp = gradient(d_dp);
        dp = abs(dp);
        [~, pos] = sort(dp);
        index = zeros(1,length(p_norm));

        k = 0;
        for i = 1:length(pos)
            if (t(pos(i)) == t(1) && sign(d_dp(1)) < 0)
                index(pos(i)) = 1;
                k = k+1;
            elseif (t(pos(i)) == t(end) && sign(d_dp(end)) > 0)
                index(pos(i)) = 1;
                k = k+1;
            elseif (sign(dd_dp(pos(i))) < 0)
                index(pos(i)) = 1;
                k = k+1;
            end

            if (k == m)
                break;
            end
        end
        
        index = logical(index);
        t_pruned = t_pruned(index);
        index = kron(index, ones(1,n));
        
    elseif (~isempty(Phi))
        t_pruned = t_pruned(logical(imp_opp));
        index = kron(imp_opp, ones(1,n));
    end

    dv = zeros(n,N);
    if (~isempty(Phi))
        A = Phi(logical(index),:);
        dv = reshape((A.'\b), n, []); 
    end

    % Output
    switch (obj.Thruster.p)
        case 'L2'
            obj.Cost = sqrt( dot(dv, dv, 1) );
        case 'L1'
            obj.Cost = sum( abs(dv), 1 );
        case 'Linfty'
            obj.Cost = max( abs(dv), [], 1 );
    end

%     k = 1;
%     for i = 1:1000:size(z,2)
%         dv = reshape(z(m+1:end,i), n, []);
%         dVs(:,k) = sqrt(dot(dv,dv,1));
%         k = k+1;
%     end
% 
%     figure 
%     stem3(1:1000:size(z,2),t,dVs, 'LineStyle','none');

    dV = zeros(n, length(t));
    for i = 1:length(t_pruned)
        dV(:, t_pruned(i) == t) = dv(:,i);
    end

    obj.Cost = dot(b,lambda);
   
    obj.Report = Output;                        % Optimization report
    obj.e(:,1) = b - M.' * reshape(dV, [], 1);  % Final missvector   
    obj.u = dV;                                 % Final rendezvous impulsive sequence
    obj.t = t;                                  % Execution times
    e = obj.e;
end