%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 11/10/25
% File: SolverADMM.m 
% Issue: 0 
% Validated: 

%% ADMM solver %%
% This function contains the main ADMM solver class

classdef SolverADMM
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
    
    % Dimensions of the problem
    properties (Access = private)
        m;                          
        n;                         
        j; 
        k;
    end

    methods
        % Constructor 
        function [obj] = SolverADMM(myObjective, myX_update, myZ_update, myRho, myA, myB, myC)
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
        function [x, z, Output] = solver(obj)
            SolverCheck = isequal( obj.A, eye(obj.m) ); 
            SolverCheck = SolverCheck && isequal( obj.B, -eye(obj.j) );
            SolverCheck = SolverCheck && isequal( obj.C, zeros(obj.m) );

            if ( SolverCheck )
                [x, z, Output] = obj.eye_solver();
            else
                [x, z, Output] = obj.general_solver();
            end
        end
    end

    methods (Access = private)
        % Solver for the simplified ADMM problem 
        [x, z, Output] = eye_solver( obj );

        % Solver for the general ADMM problem
        [x, z, Output] = general_solver( obj );
    end
end