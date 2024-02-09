% Sergio Cuevas del Valle
% Date: 25/01/23
% File: PVT_pruner.m 
% Issue: 0 
% Validated: 25/01/23

%% Primer Vector Theory Pruner %%
% This script contains the function to reduce an impulsive control law by means of Potter's generalized pruner.

% Inputs: - matrix Phi, the STM matrix of the system in time
%         - matrix B, the control matrix of the system (possibly in time)
%         - array dV, the impulsive, non-pruned, control law
%         - scalar dVmax, the maximum control authority of the thuster
%         - scalar dVmin, the minimum control authority of the thuster
%         - scalar p, the fuel-consumption proxy metric used

% Output: - array dV, the pruned impulses sequence 
%         - scalar cost, the associated p-cost of the sequence

% New versions: 

%% Functions
function [dV, cost] = PVT_pruner(Phi, B, dV, dVmax, dVmin, p, equil_flag)
    % Constants 
    m = size(Phi,1);                              % Dimension of the state space
    n = size(dV,1);                               % Dimension of the control space
    N = size(dV,2);                               % Length of the sequence

    % Sanity checks
    if (~exist('equil_flag', 'var'))
        equil_flag = false;
    end

    if (~exist('dVmax', 'var'))
        dVmax = Inf * ones(1,N);
    elseif (dVmax <= 0)
        error('Maximum control authority must be positive.')
    elseif (size(dVmax) == [1 1])
        dVmax = repmat(dVmax, 1,N);
    end

    if (~exist('dVmin', 'var'))
        dVmin = zeros(1,N);
    elseif (dVmax <= 0)
        error('Minimum control authority must be positive.')
    elseif (size(dVmin) == [1 1])
        dVmin = repmat(dVmin, 1,N);
    end

    if (~exist('p', 'var'))
        p = 'L2';
    end

    % Compute the initial cost and get the Holder conjugate norm
    switch (p)
        case 'L1'
            q = 'Linfty';
            cost = sum( abs(dV), 1 );
            num_sequence = N - (m+1);           % Number of sequence reductions

        case 'L2'
            q = 'L2';
            cost = sqrt( dot(dV, dV, 1) );
            num_sequence = N - (m+1);           % Number of sequence reductions

        case 'Linfty'
            q = 'L1';
            cost = max( abs(dV), [], 1 );
        otherwise
            error('No valid thruster configuration was selected');
    end
    cost = sum(cost);
    

    if (num_sequence >= 0)
        if (size(B,2) == size(dV,1))
            B = repmat(B, [1 N]);
        end

        % Get the initial solution
        switch (q)
            case 'Linfty'
                % Final dynamic matrix
                M = Phi(:,end-m+1:end);
            
                u = zeros(m, n * N);
                for i = 1:N
                    R = M * ( Phi(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
                    u(:,1 + n * (i-1) : n * i) = R;
                end

                A = u;

                if (any(dVmax ~= Inf) || any(dVmin > 0))       
                    A = [A zeros(m,N); ...                       % Dynamics
                         +eye(n*N) -kron(eye(N),ones(n,1));      % Epigraph form
                         -eye(n*N) -kron(eye(N),ones(n,1));      % Epigraph form
                         ]; 

                    if (any(dVmax ~= Inf))
                        A = [A; [zeros(N, n * N) eye(N)] ];
                    end

                    if (any(dVmin > 0))
                        A = [A; [zeros(N, n * N) -eye(N)] ];
                    end
                end

                qf = ones(n * N,1);

                % Equilibration
                if (equil_flag)
                    [qf, A, ~, D1, D2] = Solvers.Ruiz_equil(qf, A, 1e-6, 'R');
                    D2 = diag(D2).';
                    D1 = diag(D1).';
                else
                    D2 = ones(1, n * N);
                    D1 = ones(1, size(A,1));
                end
                
                % Initial setup
                u = A(1:m,:);
                con_vector = zeros(m,1);

                % Basis pursuit formulation
                B = reshape(dV, 1, n * N) ./ D2(1,:);

                V = zeros(2, size(B,1));
                V(1, B > 0) = +B(B > 0);
                V(2, B < 0) = -B(B < 0);

                V = [reshape(V(1,:), n, N); reshape(V(2,:), n, N)];
                Vnorm = max(abs( reshape(B, n, N) ), [], 1);

                if (any(dVmax ~= Inf))
                    
                    % Equilibration of the bounds
                    dVmax = dVmax .* D1(1, m+2*n*N+1:m+2*n*N+N);

                    if (any(Vnorm > dVmax))
                        error('Pruner cannot continue. Infeasible initial solution detected.')
                    else

                        V = [V(1:2*n,:); ...    % Solutions
                             +Vnorm; ...        % Epigraph form
                             +Vnorm; ...        % Epigraph form
                             dVmax - Vnorm];    % Slack
                    end
                end
    
                if (any(dVmin > 0))

                    % Equilibration of the bounds
                    dVmin = dVmin .* D1(1, end-N+1:end);

                    if (any(Vnorm < dVmin))
                         error('Pruner cannot continue. Infeasible initial solution detected.')
                    else
                        V = [V; ...
                             +Vnorm; ...        % Epigraph form
                             +Vnorm; ...        % Epigraph form
                             Vnorm - dVmin];    % Slack
                    end
                end
    
            case 'L2'
                Vnorm = sqrt(dot(dV, dV, 1));             % L2 norm
                
                % Final dynamic matrix
                M = Phi(:,end-m+1:end);
            
                u = zeros(m,N);
                for i = 1:N
                    R = M * ( Phi(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
                    u(:,i) = R * dV(1:n,i);
                end

                Idx = Vnorm ~= 0;                         % Non-zero elements
                u(:,Idx) = u(:,Idx) ./ Vnorm(Idx);        % Normal costs

                % Equlibration
                A = u;
                if (any(dVmax ~= Inf))       
                    A = [A; eye(N)]; 
                end

                if (any(dVmin > 0))
                    A = [A; -eye(N)];
                end

                qf = ones(N,1);
                
                % Equilibration
                if (equil_flag)
                    [qf, A, ~, D1, D2] = Solvers.Ruiz_equil(qf, A, 1e-5, false);
                    D2 = diag(D2).';
                    D1 = diag(D1).';
                else
                    D2 = ones(1,N);
                    D1 = ones(1,size(A,1));
                end

                % Initial setup
                u = A(1:m,:);
                con_vector = zeros(m,1);   % Constraint vector
                Vnorm = Vnorm ./ D2;       % Equilibration of the initial sequence
                V = Vnorm;

                if (any(dVmax ~= Inf))
                    
                    % Equilibration of the bounds
                    dVmax = dVmax .* D1(1, m+1:m+N);

                    if (any(Vnorm > dVmax))
                        error('Pruner cannot continue. Infeasible initial solution detected.')
                    else

                        V = [V(1,:); ...
                             dVmax - Vnorm];
                    end
                end
    
                if (any(dVmin > 0))

                    % Equilibration of the bounds
                    dVmin = dVmin .* D1(1, end-N+1:end);

                    if (any(Vnorm < dVmin))
                         error('Pruner cannot continue. Infeasible initial solution detected.')
                    else
                        V = [V; ...
                             Vnorm - dVmin];
                    end
                end
        
            case 'L1'
                error('Constrained Linfty rendezvous is not yet implemented');
        end

        % Prune the sequence
        for i = 1:num_sequence+1
            index = 1:m+i;
            null_flag = true;

            switch (p)
                case 'L1'                    
                    idx = logical( kron(index, ones(1,n)) );
                    
                case 'L2'
                    idx = index;

                otherwise
            end
           
            while (null_flag)
                [V(:,index), cost, null_flag] = sequence_reduction(m, n, p, q, u(:,idx), qf(idx,1), V(:,index), con_vector, dVmax(:,index), dVmin(:,index));
            end

        end
    
        % Final sequence
        switch (p)
            case 'L2'
                u = dV; 
                u(:,Idx) = u(:,Idx) ./ sqrt(dot(dV(:,Idx), dV(:,Idx), 1));
                dV = (V(1,:) .* D2) .* u; 

            case 'L1'
                dV = V(1:n,:) - V(n+1:2*n,:);
                dV = reshape(dV, 1, n*N) .* D2(1,:);
                dV = reshape(dV, n, N);

            case 'Linfty'
        end
    end

    % Clear persistent variables 
    clear PVT_pruner;
end

%% Auxiliary functions 
function [x, cost, null_flag] = sequence_reduction(m, n, p, q, u, qf, x, lambda, xmax, xmin)
    % Final indices
    switch (q)
        case 'L2'
            Indx = x(1,:) ~= 0;                 % Feasible impulses, non-vertex      
            
            % Add the constraints
            if (any(xmax ~= Inf))
                Indx = Indx & x(1,:) < xmax;
                Indx = Indx & x(2,:) ~= 0;
                xmax = xmax.';
            end

            if (any(xmin > 0))
                Indx = Indx & x(1,:) > xmin; 
                xmin = xmin.';
            end

            N = sum(Indx);                          % Number of non-zero variables
            U = u(:,Indx);                          % Considered subset of the sequence
            qf = qf(Indx,1).';                      % Cost function
            V = reshape(x(:,Indx).', 1, []);        % Get all the variables in a row
            xmax = xmax(Indx);                      % Upper bounds
            xmin = xmin(Indx);                      % Lower bounds

            % Assemble the constraint matrix 
            dim = size(x,1);
            if any(xmax ~= Inf)
                U = [U zeros(m, N); eye(N) eye(N)]; % Complete matrix
                lambda = [lambda; xmax];            % Complete independent term
                qf = [qf zeros(1,N)];               % Complete cost function
            end
    
            if any(xmin > 0)
                if (dim > 2)
                    U = [U zeros(size(U,1), N); eye(N) zeros(N) -eye(N)];
                else
                    U = [U zeros(m, N); eye(N) -eye(N)];
                end
                lambda = [lambda; xmin];
                qf = [qf zeros(1,N)];
            end

            null_flag = N > m;                      % Flag to indicate if the sequence is reducible

        case 'Linfty'
            Indx = sum(x(1:2*n,:)) ~= 0;
            N = sum(Indx);                          % Number of non-zero impulses
            null_flag = N > m;                      % Flag to indicate if the sequence is reducible

            if (null_flag)    
                % Feasible impulses, non-vertex
                U = [u -u];
                qf = [qf; qf];

                if (size(x,1) > 2 * n)
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, []); 
                        [x(2*n+1:end,:) ones(size(x,1)-2*n, 2 * n * N - size(x,2))]];
                else
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, [])];
                end
                Indx = X(1,:) ~= 0;
    
                if (xmax ~= Inf)
                    Indx = Indx & repmat( kron(x(2*n+3,:) < xmax, ones(1,n)), 1, 2);
                end
    
                if (xmin > 0)
                    Indx = Indx & repmat( kron(x(2*n+3,:) > xmin, ones(1,n)), 1, 2); 
                    Indx = Indx & repmat( kron(x(end,:) ~= 0, ones(1,n)), 1, 2);
                end
    
                U = U(:,Indx);                          % Considered subset of the sequence
                qf = qf(Indx,1).';                      % Cost function
                N = sum(Indx);
                
                % Assemble the constraint matrix 
                dim = size(x,1);
                if any(xmax ~= Inf)
                    U = [U zeros(m, 2*N); 
                        +eye(n*N) -eye(n*N) -kron([eye(N) zeros(N)],ones(n,1));
                        -eye(n*N) +eye(n*N) -kron([eye(N) zeros(N)],ones(n,1));
                         zeros(N, 2 * n * N) eye(N) eye(N)];                      % Complete matrix

                    lambda = [lambda; zeros(2 * n * N,1); xmax];            % Complete independent term
                    qf = [qf zeros(1,N)];                                   % Complete cost function
                end
    
                if any(xmin > 0)
                    if (dim > 2 * n)
                        U = [U zeros(size(U,1), N); eye(N) zeros(N) -eye(N)];
                    else
                        U = [U zeros(m, N); 
                            +eye(n*N) -eye(n*N) -kron(eye(N),ones(n,1));
                            -eye(n*N) +eye(n*N) -kron(eye(N),ones(n,1));
                            zeros(N, 2 * n * N) eye(N) -eye(N)];            % Complete matrix
                    end

                    lambda = [lambda; xmin];
                    qf = [qf zeros(1,N)];
                end

            end
            
        otherwise
    end

    if (null_flag)
       % Update of the SVD
        persistent svd_set
        persistent S
        persistent R
        persistent L

%         if isempty(svd_set)
%             [L, S, R] = svd(U);
%             svd_set = false;
%         else
%             % Rank 1 update 
%             R = [R; zeros(1,size(R,2))];
%             b = [zeros(1,size(R,2)) 1].';
%             [L1, S1, R1] = Solvers.Brand_SVD(L, S, R, U(:,end), b);
% 
%             [L, S, R] = svd(U);
%             norm(R-R1)
%             norm(L1-L)
%             norm(S1-S)
%         end

        % Compute the coordinate vector associated to a null impulse
        [L, S, R] = svd(U);
        [Lb, Ub] = updateLU(U,2);
        e = LUsolve(Lb, Ub, lambda);
        v = R(:, sum(S,1) == 0);
        v = v ./ sqrt(dot(v,v,1));
        
        c = ones(size(v,2),1);
        alpha = (e + v * c).';
        
        if (sum(lambda) == 0)
            Alpha = dot( qf, alpha );           % New cost function
            if (Alpha < 0)
                alpha = -alpha;                 % Ensure feasibility of the solution in the unconstrained case
            end
        end
    
        % Update the impulse sequence 
        beta = V ./ alpha;                      % Compute the sequence ratios
        beta_r = min(beta(alpha > 0));          % Minimum ratio

        if (beta_r > 0)
            mu = 1 - beta_r ./ beta;            % New sequence magnitudes

%             if (any(mu == 0))
%                 idc = mu == 0;
%                 for i = 1:length(mu)
%                     % Rank 1 update 
%                     if (idc(i))
%                         b =  [repmat(0, 1, i-1) 1 repmat(0, 1, size(R,2)-i)];
%                         [L, S, R] = Solvers.Brand_SVD(L, S, R, -U(:,i), b.');
%                         L = L(:,~b);
%                         S = S(~b, ~b);
%                         R = R(~b,~b);
%                     end
%                 end
%             end
            
            % New sequence and slack variables
            switch (q)
                case 'L2'
                    V = mu .* V;
                    x(:,Indx) = reshape(V, size(x(:,Indx),2), size(x,1)).';
                                       
                case 'Linfty'
                    V = mu .* V;
                    X(:,Indx) = reshape(V, size(X(:,Indx),2), size(X,1)).';

                    if (size(x,1) > 2 * n)
                    else
                        dim = size(X,2) / 2;
                        x = [reshape(X(:,1:dim), n, []); reshape(X(:,dim+1:end), n, [])];
                    end
            end
        else
            warning('Infeasibility detected. Adding coasting will not improve the solution')
        end
    end

    % Final cost and saturation
    switch (p)
        case 'L2'
            cost = x(1,:);
            null_flag = false;                      % Flag to indicate if the sequence is reducible
            
        case 'L1'
            cost = sum( abs(x(1:2*n,:)), 1 );

        case 'Linfty'
            cost = max( abs(dV(1:n,:)), [], 1 );
    end

    % Final cost
    cost = sum(cost);
end

%% Auxiliart functions 
function [x] = LUsolve(L,U,b)
    y = L\b;
    x = U\y;
end