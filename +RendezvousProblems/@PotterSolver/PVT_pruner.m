% Sergio Cuevas del Valle
% Date: 25/01/23
% File: PVT_pruner.m 
% Issue: 0 
% Validated: 25/01/23

%% Primer Vector Theory Pruner %%
% This script contains the function to reduce an impulsive control law by means of Potter's pruner

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
function [dV, cost] = PVT_pruner(Phi, B, dV, p)
    % Constants 
    m = size(Phi,1);                              % Dimension of the state space
    N = size(dV,2);                               % Length of the sequence

    % Compute the initial cost and get the Holder conjugate norm
    switch (p)
        case 'L2'
            cost = sqrt( dot(dV, dV, 1) );
            num_sequence = N - (m+1);             % Number of sequence reductions

        otherwise
            error('No valid thruster configuration was selected');
    end
    cost = sum(cost);
    
    if (num_sequence >= 0)
        if (size(B,2) == size(dV,1))
            B = repmat(B, [1 N]);
        end

        % Get the initial solution and the problem setup
        [V, ~, A, ~] = L2_preparation(Phi, B, dV);

        % Prune the sequence
        for i = 1:num_sequence+1
            idx = 1:m+i;
            [V(:,idx), cost] = RendezvousProblems.PotterSolver.sequence_reduction(m, A(:,idx), V(:,idx));
        end
    
        % Final sequence
        u = dV; 
        Idx = sqrt(dot(dV,dV,1)) ~= 0;
        u(:,Idx) = u(:,Idx) ./ sqrt(dot(dV(:,Idx), dV(:,Idx), 1));
        dV = V(1,:).* u; 
    end
end

%% Auxiliary functions 
% Prepare the L2 problem 
function [V, qf, A, b] = L2_preparation(Phi, B, dV)
    % Constants 
    m = size(Phi,1);                           % Dimension of the state space
    n = size(dV,1);                            % Dimension of the control space
    N = size(dV,2);                            % Length of the sequence

    % L2 norm
    Vnorm = sqrt(dot(dV, dV, 1)); 
    V = Vnorm;
                
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
    b = zeros(m,1);                           % Constraint vector
end