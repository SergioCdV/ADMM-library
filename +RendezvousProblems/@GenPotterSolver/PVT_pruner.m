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
    num_sequence = size(dV,2)-(size(Phi,1)+1);    % Number of sequence reductions
    sequence = dV;                                % Initial sequence

    % Sanity checks
    if (size(B,2) == size(dV,1))
        B = repmat(B, [1 size(dV,2)]);
    end

    if (~exist('dVmax', 'var'))
        dVmax = Inf;
    end

    if (~exist('dVmin', 'var'))
        dVmin = 1E-6;
    end

    if (~exist('p', 'var'))
        p = 'L2';
    end

    switch (p)
        case 'L1'
            q = 'Linfty';
        case 'L2'
            q = 'L2';
        case 'Linfty'
            q = 'L1';
        otherwise
            error('No valid thruster configuration was selected');
    end

    if (size(dV,1) ~= n)
        try 
            dV = reshape(dV, n, []);
        catch
            error('Input control dimension mismatches the initial solution');
        end
    end

    % Compute the initial cost 
    switch (p)
        case 'L2'
            cost = sqrt( dot(dV,dV,1) );
        case 'L1'
            cost = sum( abs(dV), 1 );
        case 'Linfty'
            cost = max( abs(dV), [], 1 );
    end

    cost = sum(cost);

    % Initial setup
    if (num_sequence >= 0)
        % Get the sequence 
        V = sqrt(dot(sequence, sequence,1));
        index = V > 0 & V <= dVmax;

        while (sum(index) >= m+1)
            % Reduce the sequence
            seq = sequence(:,index);
            Pos = logical( kron(index, ones(1,m)) );
            phi = Phi(:,Pos);
            Pos = logical( kron(index, ones(1,n)) );
            b = B(:,Pos);
            pos = 1:m+1:size(seq,2);

            if (length(pos) > 1)
                for i = 1:length(pos)-1
                    Pos_s = logical( kron(pos(i):pos(i+1)-1, ones(1,m)) );
                    Pos_c = logical( kron(pos(i):pos(i+1)-1, ones(1,n)) );
                    [seq(:,pos(i):pos(i+1)-1), cost] = sequence_reduction(phi(:,Pos_s) , b(:,Pos_c), seq(:,pos(i):pos(i+1)-1), dVmax, dVmin, p, q);
                end
            else
                [seq, cost] = sequence_reduction(phi, b, seq, dVmax, dVmin, p);
            end

            % Final sequence
            sequence(:,index) = seq;

            % Get the sequence 
            V = sqrt(dot(sequence, sequence,1));
            index = V > 0 & V <= dVmax;
        end
    
        % Final sequence 
        dV = sequence; 
    end
end

%% Auxiliary functions 
function [dV, cost] = sequence_reduction(Phi, B, dV, dVmax, dVmin, p, q)
    % Constants 
    m = size(B,1);      % Dimension of the state space
    n = size(dV,1);     % Dimension of the control space

    % Final dynamic matrix
    M = Phi(:,end-m+1:end);

    u = zeros(m, size(dV,2));
    for i = 1:size(dV,2)
        R = M * ( Phi(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
        u(:,i) = R * dV(:,i);
    end

    % Final globa matrix
    V = sqrt( dot(dV, dV, 1) );

    u(:, V ~= 0) = u(:, V ~= 0)./V(V ~= 0);
    index = V == 0 | V >= dVmax;
    dumb = u(:,~index);

    if (sum(~index) >= 2)
        % Compute the coordinate vector associated to a null impulse
        v = (-dumb(:,1:end-1)\dumb(:,end)).';
        dumb = [v 1];
    
        alpha = zeros(1,size(u,2));
        alpha(~index) = dumb;

        Alpha = sum(alpha);
        if (Alpha < 0)
            alpha = -alpha;                     % Ensure feasibility of the solution
        end
    
        % Update the impulse sequence 
        beta = alpha(V ~= 0) ./ V(V ~= 0);      % Compute the sequence ratios
        beta_r = max(beta);                     % Minimum ratio
        mu = 1 - beta ./ beta_r;                % New sequence magnitudes
        dV(:,V ~= 0) = mu .* dV(:,V ~= 0);      % New sequence
    end

    % Final cost
    switch (p)
        case 'L2'
            cost = sqrt( dot(dV,dV,1) );
        case 'L1'
            cost = sum( abs(dV), 1 );
        case 'Linfty'
            cost = max( abs(dV), [], 1 );
    end
    cost = sum(cost);
end