% Sergio Cuevas del Valle
% Date: 25/01/23
% File: PVT_pruner.m 
% Issue: 0 
% Validated: 25/01/23

%% Primer Vector Theory Pruner %%
% This script contains the function to reduce an impulsive control law by means of Lawden's pruner.

% Inputs: - matrix Phi, the STM matrix of the system in time
%         - matrix B, the control matrix of the system (possibly in time)
%         - array dV, the impulsive, non-pruned, control law
%         - string norm, the norm to be used

% Output: - array dV, the pruned impulses sequence 
%         - scalar cost, the associated L2 cost of the sequence

% New versions: 

%% Functions
function [dV, cost] = PVT_pruner(Phi, B, dV, norm)
    % Constants 
    m = size(B,1);                              % Dimension of the state space
    num_sequence = size(dV,2)-(size(B,1)+1);    % Number of sequence reductions

    % Initial setup
    if (num_sequence >= 0)
        sequence = dV;                      % Initial sequence

        % Get the sequence 
        V = sqrt(dot(sequence, sequence,1));
        index = V ~= 0;
        cost = 0;

        while (sum(index) >= m+1)
            % Reduce the sequence
            seq = sequence(:,index);
            phi = Phi(index,:);
            pos = 1:m+1:size(seq,2);

            if (length(pos) > 1)
                for i = 1:length(pos)-1
                    [seq(:,pos(i):pos(i+1)-1), cost] = sequence_reduction(phi(pos(i):pos(i+1)-1,:), B, seq(:,pos(i):pos(i+1)-1), norm);
                end
            else
                [seq, cost] = sequence_reduction(phi, B, seq, norm);
            end

            % Final sequence
            sequence(:,index) = seq;

            % Get the sequence 
            V = sqrt(dot(sequence, sequence,1));
            index = V ~= 0;
        end
    
        % Final sequence 
        dV = sequence; 
    else
        switch (norm)
            case 'L2'
                cost = sum(sqrt(dot(dV,dV,1)));
            case 'L1'
                cost = sum(sum(abs(dV),1));
        end
    end
end

%% Auxiliary functions 
function [dV, cost] = sequence_reduction(Phi, B, dV, norm)
    % Constants 
    m = size(B,1);                              % Dimension of the state space

    % Final matrix
    M = reshape(Phi(end,:), [m m]);

    u = zeros(m, size(dV,2));
    for i = 1:size(dV,2)
        R = M*reshape(Phi(i,:), [m m])^(-1)*B;
        u(:,i) = R*dV(:,i);
    end

    V = sqrt(dot(dV,dV,1));
    u(:, V ~= 0) = u(:, V ~= 0)./V(V ~= 0);
    index = V == 0;
    dumb = u(:,~index);

    if (sum(~index) >= 2)
        v = (-dumb(:,1:end-1)\dumb(:,end)).';
        diff = V - sqrt(3)/3 * sum(abs(dV),1);
        res = sum(v) + (1 - max([v 1]./V(V ~= 0)) * sum(diff));
        alpha_g = -sign(res);

        dumb = [v 1] * alpha_g;
    
        alpha = zeros(1,size(u,2));
        alpha(~index) = dumb;
        Alpha = sum(alpha);
    
        if (Alpha < 0)
            alpha = -alpha;
        end
    
        % New impulse sequence 
        beta = alpha(V ~= 0)./V(V ~= 0);
        beta_r = max(beta);

        check = sum(sum(abs(dV),1)) <= sqrt(3) * (sum(V)-Alpha/beta_r);

        if (~check)
            a = 1;
        end

        mu = 1./beta_r.*(beta_r-beta);
        dV(:,V ~= 0) = mu.*dV(:,V ~= 0); 
    end

    % Final cost
    switch (norm)
        case 'L2'
            cost = sum(sqrt(dot(dV,dV,1)));
        case 'L1'
            cost = sum(sum(abs(dV),1));
    end
end

