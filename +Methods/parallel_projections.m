%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/03/23
% File: parallel_projections.m 
% Issue: 0 
% Validated: 

%% Parallel projections %%
% This function contains the main PP solver class

classdef parallel_projections
    % Main properties
    properties 
        X_update;                   % Prime optimization problem solver

        % Method hyperparameters
        MaxIter = 1e3;              % Maximum number of iterations
        AbsTol = 1e-4;              % Absolute tolerance
        RelTol = 1e-2;              % Relative tolerance

        QUIET = true;               % Output results flag

        N;
        n; 
    end

    methods 
        % Constructor 
        function [obj] = parallel_projections(myX_update, myN, myn)
            obj.X_update = myX_update;
            obj.N = myN;
            obj.n = myn;
        end

        % Options
        function [obj] = optionsADMM(myMaxIter, myAbsTol, myRelTol, myQuiet)
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
        function [x,z,Output] = solver(obj)
            [x,z,Output] = general_solver(obj);
        end
    end 

    methods (Access = private)
        function [x, z, Output] = general_solver(obj)
            % Initialization 
            iter = 1; 
            GoOn = true;

            L = obj.n*obj.N;
            x = ones(L,obj.MaxIter+1);
            z = zeros(L,obj.MaxIter+1);
            u = zeros(L,1);

            tic;
            % Solver
            while (GoOn && iter <= obj.MaxIter)
                % X update
                x(:,iter+1) = feval(obj.X_update, x(:,iter), z(:,iter), u);

                % Z update 
                xm = zeros(obj.n,1);
                um = zeros(obj.n,1);
                for i = 1:obj.N
                    xm = xm + x(1+obj.n*(i-1):obj.n*i,iter+1);
                    um = um + u(1+obj.n*(i-1):obj.n*i,1);
                end
                xm = xm/obj.N;
                um = um/obj.N;

                z(:,iter+1) = repmat(xm+um,obj.N,1);

                % Compute the residuals and update the 
                u = u + (x(:,iter+1) - repmat(xm,obj.N,1));
                
                % Convergence analysis 
                Output.r_norm(iter) = norm(x(:,iter+1) - z(:,iter+1));
                Output.s_norm(iter) = norm(-(z(:,iter+1) - z(:,iter)));
            
                Output.eps_pri(iter) =  sqrt(L) * obj.AbsTol + obj.RelTol * max(norm(x(:,iter+1)));
                Output.eps_dual(iter) = sqrt(L) * obj.AbsTol + obj.RelTol * norm(u);

                if (obj.QUIET)
                    if (iter == 1)
                        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s', 'iter', 'r norm', 'eps pri', 's norm', 'eps dual');
                        fprintf('\n-----------------------------------------------------------------------------------\n');
                        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter));
                    else
                        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter));
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
            
            if (GoOn)
                x = x(:,1:iter-1); 
                z = z(:,1:iter-1);
                Output.Iterations = iter-1;
            else
                x = x(:,1:iter); 
                z = z(:,1:iter);
                Output.Iterations = iter;
            end
        end
    end
end