%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/12/23
% File: Ruiz_equilibration.m 
% Issue: 0 
% Validated: 

%% Ruiz equilibration %%
% This function contains the implementation of the modified Ruiz
% equilibration procedure for general convex cones

function [Pt, qt, At, c, D, E] = Ruiz_equilibration(P, q, A, eps)
    % Initilization 
    n = size(P);    % Dimension of the scaling matrix
    c = 1;          % Cost scaling term
    S = eye(2*n);   % Equilibration matrix 

    Pt = P;         % Equilibrated quadratic penalty matrix 
    qt = q;         % Equilibrated cost function
    At = A;         % Equilibrated equality constraint

    % Ruiz equilibration
    iter = 1; 
    maxIter = 100;
    GoOn = true;

    while (GoOn && iter < maxIter)
        % Form the M matrix 
        M = [Pt At.'; At zeros(n)]; 

        delta = 1 ./ sqrt(max(abs(M), [], 1));

        Delta = diag(delta);
        Pt = Delta * Pt; 
        qt = Delta * qt;
        At = Delta * At; 

        gamma = 1 / max([mean(max(abs(Pt), [], 1)), norm(qt, 'inf')]);
        Pt = gamma * Pt;
        qt = gamma * qt;
        S = Delta * S;
        c =  gamma *  c;

        % Convergence analysis 
        if (norm(1-delta, 'inf') < eps)
            GoOn = false;
        else
            iter = iter + 1;
        end
    end

    % Final scaling matrices
    D = S(1:n,1:n);
    E = S(n+1:end, n+1:end);
end