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
    n = size(P,1);    % Dimension of the scaling matrix
    N = size(A,1);    % Dimension of the equality matrix
    S = eye(n + N);   % Equilibration matrix 
    c = 1;            % Cost scaling term

    Pt = P;           % Equilibrated quadratic penalty matrix 
    qt = q;           % Equilibrated cost function
    At = A;           % Equilibrated equality constraint

    % Ruiz equilibration
    iter = 1; 
    maxIter = 100;
    GoOn = true;

    while (GoOn && iter < maxIter)
        % Form the M matrix 
        M1 = [Pt; At];

        delta = [max(abs(M1), [], 1) max(abs(At), [], 2).'];
        delta(delta < 1e-5) = ones(1, length(delta(delta < 1e-5)));
        Delta = diag( 1 ./ sqrt(delta) );

        S = Delta * S;
        D = Delta(1:n, 1:n);
        E = Delta(n+1:end, n+1:end);

        % Ruiz equilibration
        Pt = D * (Pt * D); 
        At = E * (At * D); 
        qt = D * qt;

        % Cost equilibration step
        gamma = max([mean(max(abs(Pt), [], 1)), norm(qt, 'inf')]);
        if (gamma < 1e-5)
            gamma = 1;
        else
            gamma = 1 / gamma;
        end

        Pt = gamma * Pt;
        qt = gamma * qt;
        c =  gamma * c;

        % Convergence analysis 
        if (norm(1-delta, 'inf') < eps)
            GoOn = false;
        end
        
        iter = iter + 1;
    end

    % Final scaling matrices
    D = S(1:n, 1:n);
    E = S(n+1:end, n+1:end);
end