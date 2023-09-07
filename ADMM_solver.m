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

        x;                          % First consensus variable
        z;                          % Second consensus variable
        u;                          % Lagrange penalizer

        % Method hyperparameters
        alpha = 0;                  % Relaxation coefficient
        MaxIter = 1e4;              % Maximum number of iterations
        AbsTol = 1e-9;              % Absolute tolerance
        RelTol = 1e-6;              % Relative tolerance

        QUIET = true;               % Output results flag
    end

    properties (Access = private)
        m;
        n; 
        j; 
        k;
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
            obj.n = size(obj.A,2);
            obj.j = size(obj.B,1);
            obj.k = size(obj.B,2);

            if (obj.j ~= obj.m)
                error('No valid matrix dimensions were introduced.');
            else
                obj.x = zeros(obj.n,obj.MaxIter+1);
                obj.z = zeros(obj.k,obj.MaxIter+1);
                obj.u = zeros(obj.m,1);
            end
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
        function [x,z,Output] = solver(obj)
            if ( isequal(obj.A,eye(size(obj.A,1))) && isequal(obj.B,-eye(size(obj.B,1))) && isequal(obj.C, zeros(size(obj.C,1),1)) )
                solver_type = 'Eye';
            else
                solver_type = 'General';
            end
 
            switch (solver_type)
                case 'Eye'
                    [x,z,Output] = eye_solver(obj);
                otherwise
                    [x,z,Output] = general_solver(obj);
            end
        end
    end

    methods (Access = private)
        function [x, z, Output] = general_solver(obj)
            % Initialization 
            iter = 1; 
            GoOn = true;            
            tic;

            % Solver
            while (GoOn && iter <= obj.MaxIter)
                % X update
                obj.x(:,iter+1) = feval(obj.X_update, obj.x(:,iter), obj.z(:,iter), obj.u);

                % Z update with relaxation
                zold = obj.B * obj.z(:,iter) - obj.C;
                xh = obj.alpha * obj.A * obj.x(:,iter+1) - (1-obj.alpha) * zold;
                obj.z(:,iter+1) = feval(obj.Z_update, xh, obj.z(:,iter), obj.u);

                % Compute the residuals and update the 
                obj.u = obj.u + (xh + obj.B * obj.z(:,iter+1) - obj.C);
                
                % Convergence analysis 
                Output.objval(iter) = feval(obj.objective, obj.x(:,iter+1), obj.z(:,iter+1));
                Output.r_norm(iter) = norm(obj.A * obj.x(:,iter+1) + obj.B * obj.z(:,iter+1) - obj.C);
                Output.s_norm(iter) = norm(obj.rho * obj.A.' * obj.B * (obj.z(:,iter+1) - obj.z(:,iter)));
            
                Output.eps_pri(iter) =  sqrt(obj.n) * obj.AbsTol + obj.RelTol * max([norm(obj.C) norm(obj.A * obj.x(:,iter+1)), norm(obj.B * obj.z(:,iter+1))]);
                Output.eps_dual(iter) = sqrt(obj.m) * obj.AbsTol + obj.RelTol * norm(obj.rho * obj.A.' * obj.u);

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
            
            if (GoOn)
                x = obj.x(:,1:iter-1); 
                z = obj.z(:,1:iter-1);
                Output.Iterations = iter-1;
            else
                x = obj.x(:,1:iter); 
                z = obj.z(:,1:iter);
                Output.Iterations = iter;
            end
        end

        function [x, z, Output] = eye_solver(obj)
            % Initialization 
            iter = 1; 
            GoOn = true;

            tic;
            % Solver
            while (GoOn && iter <= obj.MaxIter)
                % X update
                obj.x(:,iter+1) = feval(obj.X_update, obj.x(:,iter), obj.z(:,iter), obj.u);

                % Z update with relaxation
                zold = -obj.z(:,iter);
                xh = obj.alpha * obj.x(:,iter+1) - (1-obj.alpha) * zold;
                obj.z(:,iter+1) = feval(obj.Z_update, xh, obj.z(:,iter), obj.u);

                % Compute the residuals and update the 
                obj.u = obj.u + (xh - obj.z(:,iter+1));
                
                % Convergence analysis 
                Output.objval(iter) = feval(obj.objective, obj.x(:,iter+1), obj.z(:,iter+1));
                Output.r_norm(iter) = norm(obj.x(:,iter+1) - obj.z(:,iter+1));
                Output.s_norm(iter) = norm(-obj.rho * (obj.z(:,iter+1) - obj.z(:,iter)));
            
                Output.eps_pri(iter) =  sqrt(obj.n) * obj.AbsTol + obj.RelTol * max([norm(obj.x(:,iter+1)), norm(obj.z(:,iter+1))]);
                Output.eps_dual(iter) = sqrt(obj.m) * obj.AbsTol + obj.RelTol * norm(obj.rho * obj.u);

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
            
            if (GoOn)
                x = obj.x(:,1:iter-1); 
                z = obj.z(:,1:iter-1);
                Output.Iterations = iter-1;
            else
                x = obj.x(:,1:iter); 
                z = obj.z(:,1:iter);
                Output.Iterations = iter;
            end
        end
    end
end