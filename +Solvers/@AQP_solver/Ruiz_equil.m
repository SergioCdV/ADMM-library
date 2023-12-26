%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/12/23
% File: Ruiz_equil.m 
% Issue: 0 
% Validated: 

%% Ruiz equilibration %%
% This function contains the implementation of the modified Ruiz
% equilibration procedure for general convex cones

function [qt, At, c, D2, D1] = Ruiz_equil(q, A, eps)
    % Initilization 
    N = size(A,1);    % Dimension of the equality matrix
    n = size(A,2);    % Dimension of the equality matrix 
    D1 = eye(N);      % Equilibration matrix 
    D2 = eye(n);      % Equilibration matrix
    c = 1;            % Cost scaling term

    qt = q;           % Equilibrated cost function
    At = A;           % Equilibrated equality constraint

    % Ruiz equilibration
    iter = 1; 
    maxIter = 100;
    GoOn = true;

    while (GoOn && iter < maxIter)
        delta = [max(abs(At), [], 2).' max(abs(At), [], 1)];
        delta(delta < 1e-5) = ones(1, length(delta(delta < 1e-5)));
        Delta = diag( 1 ./ sqrt(delta) );

%         E = Delta(1:N, 1:N);
        E = D1;
        D = Delta(N+1:end, N+1:end);

        % Ruiz equilibration
        At = E * (At * D); 
        qt = D * qt;
%         D1 = D1 * E;
        D2 = D2 * D;

        % Cost equilibration step
        gamma = norm(qt, 'inf');
        if (gamma < 1e-5)
            gamma = 1;
        else
            gamma = 1 / gamma;
        end

        qt = gamma * qt;
        c =  gamma * c;

        % Convergence analysis 
        if (norm(1-delta, 'inf') < eps)
            GoOn = false;
        end
        
        iter = iter + 1;
    end
end