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
function [dV, cost] = PVT_pruner(Phi, B, dV, dVmax, dVmin, p)
    % Constants 
    m = size(Phi,1);                              % Dimension of the state space
    n = size(dV,1);                               % Dimension of the control space
    N = size(dV,2);                               % Length of the sequence
    num_sequence = size(dV,2)-(size(Phi,1)+1);    % Number of sequence reductions

    % Sanity checks
    if (num_sequence >= 0)
        if (size(B,2) == size(dV,1))
            B = repmat(B, [1 size(dV,2)]);
        end
    
        if (~exist('dVmax', 'var'))
            dVmax = Inf;
        elseif (dVmax <= 0)
            error('Maximum control authority must be positive.')
        end
    
        if (~exist('dVmin', 'var'))
            dVmin = 0;
        elseif (dVmax <= 0)
            error('Minimum control authority must be positive.')
        end
    
        if (~exist('p', 'var'))
            p = 'L2';
        end
    
        % Compute the initial cost and get the Holder conjugate norm
        switch (p)
            case 'L1'
                q = 'Linfty';
                cost = sum( abs(dV), 1 );
            case 'L2'
                q = 'L2';
                cost = sqrt( dot(dV, dV, 1) );
            case 'Linfty'
                q = 'L1';
                cost = max( abs(dV), [], 1 );
            otherwise
                error('No valid thruster configuration was selected');
        end
        cost = sum(cost);
    
        % Get the initial solution
        switch (q)
            case 'L1'
                V = max(abs(dV), [], 1);
    
            case 'L2'
                V = sqrt(dot(dV, dV, 1));                   % L2 norm
                
                % Final dynamic matrix
                M = Phi(:,end-m+1:end);
            
                u = zeros(m,N);
                for i = 1:N
                    R = M * ( Phi(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
                    u(:,i) = R * dV(1:n,i);
                end

                Idx = V ~= 0;                         % Non-zero elements
                u(:,Idx) = u(:,Idx) ./ V(Idx);        % Normal costs

                % Initial setup
                con_vector = zeros(m,1);
                
                if (dVmax ~= Inf)
                    if (any(V > dVmax))
                        error('Infeasible initial solution')
                    else
                        V = [V; dVmax - V];
                    end
                end
    
                if (dVmin > 0)
                    if (any(V < dVmin))
                        error('Infeasible initial solution')
                    else
                        V = [V; V - dVmin];
                    end
                end
        
            case 'Linfty'
                error('Constrained Linfty rendezvous is not yet implemented');
                V = max(abs(sequence(1:n,:)), [], 1);
        end

        % Prune the sequence 
        for i = 1:num_sequence+1
            index = 1:m+i;
            [V(:,index), cost] = sequence_reduction(m, n, N, p, q, u(:,index), V(:,index), con_vector, dVmax, dVmin);
        end
    
        % Final sequence
        switch (q)
            case 'L2'
                u = dV; 
                u(:,Idx) = u(:,Idx) ./ sqrt(dot(dV(:,Idx), dV(:,Idx), 1));
                dV = V(1,:) .* u; 
            case 'L1'
            case 'Linfty'
        end
    end
end

%% Auxiliary functions 
function [x, cost] = sequence_reduction(m, n, N, p, q, u, x, lambda, xmax, xmin)
    % Final indices
    switch (q)
        case 'L2'
            Indx = x(1,:) ~= 0;               % Non-zero impulses
            N = sum(Indx);                    % Number of non-zero variables
            U = u(:,Indx);                    % Considered subset of the sequence
            V = reshape(x(:,Indx).', 1, []);  % Get all the variables in a column

            if (size(x,1) == 2)
                U = [U zeros(m, N); eye(N) eye(N)];
                lambda = [lambda; xmax * ones(N,1)];
            end

            if (size(x,2) == 3)
                U = [U zeros(m, N); eye(N) -eye(N)];
                lambda = [lambda; xmin * ones(N,1)];
            end

        otherwise
    end

    if (size(U,2) > m)
        % Compute the coordinate vector associated to a null impulse
        lmax = sum(lambda);
        e = lambda - U(:,end) * lmax;
        alpha = ( U(:,1:end-1) \ e ).';
        alpha = [alpha lmax];

        Alpha = sum(alpha(1:size(U,2)));
        if (Alpha < 0)
            alpha = -alpha;                     % Ensure feasibility of the solution
        end
    
        % Update the impulse sequence 
        beta = alpha ./ V;                      % Compute the sequence ratios
        beta_r = max(beta(alpha > 0));          % Minimum ratio

        if (beta_r > 0)
            mu = 1 - beta ./ beta_r;            % New sequence magnitudes
            
            % New sequence and slack variables
            switch (q)
                case 'L2'
                    V = mu .* V;
                    x(:,Indx) = reshape(V, size(x(:,Indx),2), size(x,1)).';
                    
                    idc = sum(x,1) > xmax;
                    x(:,idc) = xmax * x(:,idc) ./ sum(x(:,idc), 1);
                                       
                case 'Linfty'
            end
        else
            warning('Infeasibility detected. Adding coasting does not improve the solution')
        end
    end

    % Final cost and saturation
    switch (p)
        case 'L2'
            cost = x(1,:);
            
        case 'L1'
            cost = sum( abs(dV(1:n,:)), 1 );

        case 'Linfty'
            cost = max( abs(dV(1:n,:)), [], 1 );
    end
    a = max(cost + x(end,:))
    if (a > xmax)
        a = 1;
    end

    % Final cost
    cost = sum(cost);
end