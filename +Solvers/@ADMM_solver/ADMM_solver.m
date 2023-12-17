%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 16/12/23
% File: ADMM_solver.m 
% Issue: 0 
% Validated: 

%% ADMM solver %%
% This function contains the main general ADMM solver class

classdef ADMM_solver
    % Main properties
    properties 
        % User-defined 
        A;                          % Linear matrix of the constraint between x and z
        B;                          % Linear matrix of the constraint between x and z
        C;                          % Linear vector of the constraint between x and z

        rho = 0.1;                  % Augmented Lagrangian parameter 

        objective;                  % Objective function
        X_update;                   % Primary optimization problem solver
        Z_update;                   % Secondary optimization problem solver

        x;                          % First consensus variable
        z;                          % Second consensus variable
        u;                          % Lagrange penalizer

        % Method hyperparameters
        alpha = 1.6;                % Relaxation coefficient

        MaxIter = 1e4;              % Maximum number of iterations
        AbsTol = 1e-9;              % Absolute tolerance
        RelTol = 1e-6;              % Relative tolerance

        QUIET = true;               % Output results flag

        % Scaling (matrix equilibration)
        D;                          % Scaling sub-matrix
        E;                          % Scaling sub-matrix

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
                solver_type = 0;
            else
                solver_type = 1;
            end
            
            [x, z, Output] = solve(obj, solver_type);
        end
    end

    methods (Access = private)
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

        % General and eye solvers
        [x, z, Output] = solve(obj, solver_type);  
    end

    methods (Static, Access = private)
        [D, E] = Ruiz_equilibration(obj);   % Matrix equilibration
    end
end