%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 27/01/24
% File: sequence_reduction.m 
% Issue: 0 
% Validated: 

%% Sequence reduction %% 
% Implementation of the main Generalized Potter routine to concentrate
% candidate impulsive control plans 

% Inputs:  - m, scalar, the state space dimension 
%          - n, scalar, the control space dimension
%          - p, scalar, the p-norm of the control vector to be used in the
%            cost function (1, 2, infty)
%          - q, scalar, the q-norm of the control vector to be used in the
%            constraints (1, 2, infty)
%          - u, array, the problem's linear matrix to be used to solve the
%            problem
%          - qf, array/vector, the linear cost function of the problem
%          - x, vector, the linear solution to be concentrated

% Outputs: - x, the concentrated linear solution 
%          - cost, the lq norm of the concentrated solution 
%          - null_flag, a boolean variable to acknowledge if the solutions
%            can be further concentrated

function [x, cost] = sequence_reduction(m, u, x)
    % Final indices
    Indx = x(1,:) ~= 0;                     % Feasible impulses, non-vertex      
    N = sum(Indx);                          % Number of non-zero variables
    V = reshape(x(:,Indx).', 1, []);        % Get all the variables in a row
    null_flag = N > m;                      % Flag to indicate if the sequence is reducible

    % General selection of equations
    U = u(:,Indx);                          % Considered subset of the sequence

    if (null_flag)
        % Compute the coordinate vector associated to a null impulse
        b = -U(:,end);
        alpha = [ (U(:,1:end-1)\b).' 1 ];
        
        Alpha = sum( alpha, 2 );            % New cost function
        if (Alpha < 0)
            alpha = -alpha;                 % Ensure feasibility of the solution in the unconstrained case
        end
    
        % Update the impulse sequence 
        beta = alpha ./ V;                  % Compute the sequence ratios
        beta_r = max(beta);                 % Minimum ratio
        mu = (1 - beta ./ beta_r);          % New sequence magnitudes

        V = mu .* V;
        x(:,Indx) = reshape(V, size(x(:,Indx),2), size(x,1)).';
    end

    % Final cost and saturation
    cost = sum( x(1,:) );
end