%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/12/23
% File: Ruiz_equil.m 
% Issue: 0 
% Validated: 

%% Ruiz equilibration %%
% This function contains the implementation of the modified Ruiz
% equilibration procedure for general convex cones

function [qt, At, c, D1, D2] = Ruiz_equil(q, A, eps, lr_flag)
    % Initilization 
    N = size(A,1);    % Dimension of the equality matrix
    n = size(A,2);    % Dimension of the equality matrix 
    D1 = ones(1,N);   % Equilibration matrix 
    D2 = ones(1,n);   % Equilibration matrix
    c = 1;            % Cost scaling term

    qt = q;           % Equilibrated cost function
    At = A;           % Equilibrated equality constraint

    if (~exist('lr_flag', 'var'))
        lr_flag = false;
    end

    % Ruiz equilibration
    iter = 1; 
    maxIter = 100;
    GoOn = true;

    while (GoOn && iter < maxIter)
        if (~lr_flag)
            delta = [max(abs(At), [], 2).' max(abs(At), [], 1)];
            delta(delta == 0) = ones(1, length(delta(delta == 0)));
            Delta = 1 ./ sqrt(delta);

            E = Delta(1,1:N);
            D = Delta(1,N+1:end);
            D1 = D1 .* E;
            D2 = D2 .* D;
        else
            switch (lr_flag)
                case 'L'
                    delta = max(abs(At), [], 2).';
                    delta(delta == 0) = ones(1, length(delta(delta == 0)));
                    E = 1 ./ sqrt(delta);
                    D = D2;
                    D1 = D1 .* E;
                case 'R'
                    delta = max(abs(At), [], 1);
                    delta(delta == 0) = ones(1, length(delta(delta == 0)));
                    D = 1 ./ sqrt(delta);
                    E = D1;
                    D2 = D2 .* D;
                otherwise
                    error('No valid equilibration was selected. Abort')
            end
        end

        % Ruiz equilibration
        At = diag(E) * At * diag(D); 
        qt = diag(D) * qt;

        % Cost equilibration step
        gamma = norm(qt, 'inf');
        if (gamma == 0)
            gamma = 1;
        else
            gamma = 1 / gamma;
        end

        qt = gamma * qt;
        c =  gamma * c;

        % Convergence analysis 
        if (max( abs(1-delta) ) < eps)
            GoOn = false;
        end
        
        iter = iter + 1;
    end
end