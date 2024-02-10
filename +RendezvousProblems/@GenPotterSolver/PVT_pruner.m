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

    if (dVmax < dVmin)
        error('The minimum control authority is greater than its upper bound.')
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

        % Get the initial solution and the problem setup
        switch (q)
            case 'Linfty'
                % L1 problem
                [V, qf, A, b, D2] = L1_preparation(Phi, B, dV, dVmax, dVmin, equil_flag);
    
            case 'L2'
                % L2 problem
                [V, qf, A, b, D2] = L2_preparation(Phi, B, dV, dVmax, dVmin, 1);
        
            case 'L1'
                error('Constrained Linfty rendezvous is not yet implemented');
        end

        % Prune the sequence
        for i = 1:num_sequence+1
            index = 1:m+i;
            null_flag = true;

            switch (p)
                case 'L1'        
                    idx = zeros(1,length(index) * n);
                    for j = 1:length(index)
                        idx(1 + n*(j-1):n*j) = n * [index(j) index(j) index(j)] - (n-1:-1:0);
                    end

                    if     (size(V,1) == 4*n + 3)
                        idx = [idx, idx + n * N, index + 2 * n * N, idx + 2 * n * N + N, idx + 3 * n * N + N, index + 4 * n * N + N, index + 4 * n * N + 2 * N];

                    elseif (size(V,1) == 4*n + 2)
                        idx = [idx, idx + n * N, index + 2 * n * N, idx + 2 * n * N + N, idx + 3 * n * N + N, index + 4 * n * N + N];

                    else
                        idx = [idx, idx + n * N];
                    end
                    
                case 'L2'
                    if (size(V,1) == 3)
                        idx = [index, index + N, index + 2 * N];
                    elseif (size(V,1) == 2)
                        idx = [index, index + N];
                    else
                        idx = index;
                    end

                otherwise
            end
           
            while (null_flag)
                [V(:,index), cost, null_flag] = RendezvousProblems.GenPotterSolver.sequence_reduction(m, n, p, q, A(:,idx), qf(idx,1), V(:,index), b, ...
                                                                                                      dVmax(:,index), dVmin(:,index));
            end

        end
    
        % Final sequence
        switch (p)
            case 'L2'
                u = dV; 
                Idx = sqrt(dot(dV,dV,1)) ~= 0;
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
% Prepare the L2 problem 
function [V, qf, A, b, D2] = L2_preparation(Phi, B, dV, dVmax, dVmin, equil_flag)
    % Constants 
    m = size(Phi,1);                              % Dimension of the state space
    n = size(dV,1);                               % Dimension of the control space
    N = size(dV,2);                               % Length of the sequence

    % L2 norm
    Vnorm = sqrt(dot(dV, dV, 1));             
                
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
    qf = ones(N,1);                           % Problem's cost function
    A = u;                                    % Constrains matrix
    if (any(dVmax ~= Inf))                    % Maximum control authority constrain           
        A = [A; eye(N)]; 
    end

    if (any(dVmin > 0))                       % Minimum control authority constrain                     
        A = [A; -eye(N)];
    end
    
    % Equilibration
    if (equil_flag)
        % Ruiz equillibration
        [qf, A, ~, D1, D2] = Solvers.Ruiz_equil(qf, A, 1e-5, false);
        D2 = diag(D2).';
        D1 = diag(D1).';
    else
        % No equilibration
        D2 = ones(1,N);
        D1 = ones(1,size(A,1));
    end

    % Initial setup
    A = A(1:m,:);
    b = zeros(m,1);            % Constraint vector
    Vnorm = Vnorm ./ D2;       % Equilibration of the initial sequence
    V = Vnorm;

    if (any(dVmax ~= Inf))
        
        % Equilibration of the bounds
        dVmax = dVmax .* D1(1, m+1:m+N);

        if (any(Vnorm > dVmax))
            error('Pruner cannot continue. Infeasible initial solution detected.')
        else
            
            % Augmented decision vector: states + slacks
            V = [V(1,:); ...
                 dVmax - Vnorm];
            
            % Additional constraining
            A = [A zeros(m, N); eye(N) eye(N)];     % Complete matrix
            b = [b; dVmax.'];                       % Complete independent term
            qf = [qf; zeros(N,1)];                  % Complete cost function
        end
    end
    
    if (any(dVmin > 0))

        % Equilibration of the bounds
        dVmin = dVmin .* D1(1, end-N+1:end);

        % Augmented decision vector: states + slacks
        if (any(Vnorm < dVmin))
             error('Pruner cannot continue. Infeasible initial solution detected.')
        else
            V = [V; ...
                 Vnorm - dVmin];

            % Additional constraining
            if (size(V,1) == 3)
                A = [A zeros(size(A,1), N); eye(N) zeros(N) -eye(N)];
            else
                A = [A zeros(m, N); eye(N) -eye(N)];
            end

            b = [b; dVmin.'];                       % Complete independent term
            qf = [qf; zeros(N,1)];                  % Complete cost function
        end
    end
end

% Prepare the L1 problem
function [V, qf, A, b, D2] = L1_preparation(Phi, B, dV, dVmax, dVmin, equil_flag)
    % Constants 
    m = size(Phi,1);                              % Dimension of the state space
    n = size(dV,1);                               % Dimension of the control space
    N = size(dV,2);                               % Length of the sequence

    % Final dynamic matrix
    M = Phi(:,end-m+1:end);

    u = zeros(m, n * N);
    for i = 1:N
        R = M * ( Phi(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
        u(:,1 + n * (i-1) : n * i) = R;
    end

    A = u;                                           % Constrains matrix
    qf = ones(n * N,1);                              % Linear cost function

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
    b = zeros(m,1);     % Independent term

    % Basis pursuit formulation
    B = reshape(dV, 1, n * N) ./ D2(1,:);

    V = zeros(2, size(B,1));
    V(1, B > 0) = +B(B > 0);
    V(2, B < 0) = -B(B < 0);

    V = [reshape(V(1,:), n, N); reshape(V(2,:), n, N)];     % Basis pursuit formulation
    A = [A(1:m,1:n*N) -A(1:m,1:n*N)];                       % Basis pursuit formulation of the constraints
    qf = [qf; qf];                                          % Linear cost function

    Vnorm = max(abs( reshape(B, n, N) ), [], 1);

    if (any(dVmax ~= Inf))
        
        % Equilibration of the bounds
        dVmax = dVmax .* D1(1, m+2*n*N+1:m+2*n*N+N);

        if (any(Vnorm > dVmax))
            error('Pruner cannot continue. Infeasible initial solution detected.')
        else

            V = [V; ...
                 +Vnorm; ...                        % Epigraph form
                 +0*Vnorm - 0*dV; ...                   % Epigraph slacks
                 +0*Vnorm + 0*dV; ...                   % Epigraph slacks
                 +dVmax - Vnorm];                   % Slack
            
            % Assemble the constraint matrix
            A = [A zeros(m, N+2*n*N+N); 
                 +eye(n*N) -eye(n*N) kron(-eye(N),ones(n,1)) +eye(n*N) zeros(n*N) kron(zeros(N),ones(n,1));
                 -eye(n*N) +eye(n*N) kron(-eye(N),ones(n,1)) zeros(n*N) +eye(n*N) kron(zeros(N),ones(n,1));
                  zeros(N, 2 * n * N) eye(N) zeros(N,2*n*N) eye(N)];         

            b = [b; zeros(2 * n * N,1); dVmax.'];                           % Complete independent term
            qf = [qf; zeros(N+2*n*N+N,1)];                                  % Complete cost function
        end
    end

    if (any(dVmin > 0))

        % Equilibration of the bounds
        dVmin = dVmin .* D1(1, end-N+1:end);

        if (any(Vnorm < dVmin))
             error('Pruner cannot continue. Infeasible initial solution detected.')
        else
            % Assemble the constraint matrix 
            dim = size(V,1);
        
            % Complete matrix
            if (dim > 2 * n)
                A = [A zeros(size(A,1), N); 
                     zeros(N, 2 * n * N) eye(N) zeros(N,2*n*N + N) -eye(N)];
                
                V = [V; ...
                    Vnorm - dVmin];    % Slack
            else
                A = [A zeros(m, N+2*n*N+N); 
                     +eye(n*N) -eye(n*N) kron(-eye(N),ones(n,1)) +eye(n*N) zeros(n*N) kron(zeros(N),ones(n,1));
                     -eye(n*N) +eye(n*N) kron(-eye(N),ones(n,1)) zeros(n*N) +eye(n*N) kron(zeros(N),ones(n,1));
                      zeros(N, 2 * n * N) eye(N) zeros(N,2*n*N) -eye(N)];     

                V = [V; ...
                     +Vnorm; ...                        % Epigraph form
                     +Vnorm - dV; ...                   % Epigraph slacks
                     +Vnorm + dV; ...                   % Epigraph slacks
                     +Vnorm - dVmin];                   % Slack

                b = [b; zeros(2 * n * N,1)];            % Augment the independent term in case it has not been done yet
                qf = [qf; zeros(2*n*N+N,1)];            % Complete cost function
            end

            b = [b; dVmin.'];                           % Complete independent term
            qf = [qf; zeros(N,1)];                      % Complete cost function
        end
    end
end