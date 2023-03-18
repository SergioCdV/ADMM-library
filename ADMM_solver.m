%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% ADMM solver %%
% This function contains the main ADMM solver class

classdef ADMM_solver
    % Main properties
    properties 
        % User-defined 
        A;                          % Linear matrix of the constraint between x and z
        B;                          % Linear matrix of the constraint between x and z
        C;                          % Linear term of the constraint between x and z

        rho;                        % Augmented Lagrangian parameter 

        objective;                  % Objective function
        X_update;                   % Prime optimization problem solver
        Z_update;                   % Second optimization problem solver

        % Method hyperparameters
        alpha = 0;                  % Relaxation coefficient
        MaxIter = 1e3;              % Maximum number of iterations
        AbsTol = 1e-4;              % Absolute tolerance
        RelTol = 1e-2;              % Relative tolerance

        QUIET = true;               % Output results flag
    end

    properties (Access = private)
        m;
        n; 
    end

    methods
        % Constructor 
        function [obj] = ADMM_solver(myObjective, myX_update, myZ_update, myRho, myA, myB, myC)
            obj.rho = myRho;
            obj.A = myA;
            obj.B = myB;
            obj.C = myC;

            obj.objective = myObjective;
            obj.X_update = myX_update; 
            obj.Z_update = myZ_update;

            % Initialization
            obj = obj.initADMM();
        end

        % Initialization 
        function [obj] = initADMM(obj)
            % Dimensions of the problem
            obj.m = size(obj.A,1);
            obj.n = size(obj.B,1);
        end

        % Options
        function [obj] = optionsADMM(myAlpha, myMaxIter, myAbsTol, myRelTol, myQuiet)
            if (~isempty(myAlpha))
                obj.alpha = myAlpha;
            end

            if (~isempty(myMaxIter))
                obj.MaxIter = myMaxIter;
            end

            if (~isempty(myAbsTol))
                obj.AbsTol = myAbsTol;
            end

            if (~isempty(myRelTol))
                obj.RelTol = myRelTol;
            end

            if (~isempty(myQuiet))
                obj.QUIET = myQuiet;
            end
        end

        % Main solver
        function [x, z, Output] = solver(obj)
            % Initialization 
            iter = 1; 
            GoOn = true;

            x = zeros(obj.n,obj.MaxIter+1);
            z = zeros(obj.n,obj.MaxIter+1);
            u = zeros(obj.n,1);
            
            tic;
            % Solver
            while (GoOn && iter <= obj.MaxIter)
                % X update
                x(:,iter+1) = feval(obj.X_update, x(:,iter), z(:,iter), u);

                % Z update with relaxation
                zold = obj.B * z(:,iter) - obj.C;
                xh = obj.alpha * obj.A * x(:,iter+1) - (1-obj.alpha) * zold;
                z(:,iter+1) = feval(obj.Z_update, xh, z(:,iter), u);

                % Compute the residuals and update the 
                u = u + (xh - z(:,iter+1));
                
                % Convergence analysis 
                Output.objval(iter) = feval(obj.objective, x(:,iter+1), z(:,iter+1));
                Output.r_norm(iter) = norm(obj.A * x(:,iter+1) + obj.B * z(:,iter+1) - obj.C);
                Output.s_norm(iter) = norm(obj.rho * obj.A.' * obj.B * (z(:,iter+1) - z(:,iter)));
            
                Output.eps_pri(iter) =  sqrt(obj.n) * obj.AbsTol + obj.RelTol * max([norm(obj.C) norm(obj.A * x(:,iter+1)), norm(obj.B * z(:,iter+1))]);
                Output.eps_dual(iter) = sqrt(obj.n) * obj.AbsTol + obj.RelTol * norm(obj.rho * obj.A.' * u);

                if (obj.QUIET)
                    if (iter == 1)
                        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', 'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
                        fprintf('-----------------------------------------------------------------------------------\n');
                        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter), Output.objval(iter));
                    else
                        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter), Output.objval(iter));
                    end
                end

                % Update the state machine
                if (Output.r_norm(iter) < Output.eps_pri(iter) && Output.s_norm(iter) < Output.eps_dual(iter))
                    GoOn = false;
                else
                    iter = iter + 1;
                end
            end

            % Final results 
            Output.Time = toc;
            Output.Result = ~GoOn; 
            Output.Iterations = iter;
            
            if (GoOn)
                x = x(:,1:iter-1); 
                z = z(:,1:iter-1);
            else
                x = x(:,1:iter); 
                z = z(:,1:iter);
            end
        end
    end
end