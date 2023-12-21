%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 16/12/23
% File: AQP_solver.m 
% Issue: 0 
% Validated: 

%% AQP solver %%
% This function contains a specific ADMM solver for QP problems

classdef AQP_solver
    % Main properties
    properties 
        % User-defined 
        ConeProj;                   % Secondary optimization problem solver (projection onto the C cone)
        CheckCone;                  % Function to check if a given vector belongs to the cone

        x;                          % First consensus variable
        z;                          % Second consensus variable
        u;                          % Lagrange penalizer

        A;                          % Equality constraint matrix
        P;                          % Quadratic cost function
        q;                          % Linear cost function

        At;                         % Scaled equality constraint matrix
        Pt;                         % Scaled quadratic cost function
        qt;                         % Scaled linear function

        % Method hyperparameters
        rho = 0.1;                  % Augmented Lagrangian parameter 
        alpha = 1.6;                % Relaxation coefficient

        MaxIter = 1e4;              % Maximum number of iterations
        AbsTol = 1e-9;              % Absolute tolerance
        RelTol = 1e-6;              % Relative tolerance

        QUIET = true;               % Output results flag

        % Scaling (matrix equilibration)
        eps_equil = 1e-3;           % Equilibration tolerance
        c;                          % Cost function scaling
        D;                          % Scaling sub-matrix
        E;                          % Scaling sub-matrix

    end

    properties (Access = private)
        m;
        n;
    end

    methods
        % Constructor 
        function [obj] = AQP_solver(myConeCheck, myConeProjec, myP, myQ, myA)
            % Function handles
            obj.ConeProj = myConeProjec;
            obj.CheckCone = myConeCheck;

            % Relevant matrices and vectors
            obj.A = myA;        % Equality matrix
            obj.P = myP;        % Quadratic penalty
            obj.q = myQ;        % Cost function

            % Initialization
            obj = obj.initAQP();
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
            % Solve the problem
            [x, z, Output] = obj.quadratic_solve();

            % De-equilibration 
            x = obj.D * x;
            z = obj.D * z;
            Output.objval = Output.objval / obj.c;
        end
    end

    methods (Access = private)
        % Initialization 
        function [obj] = initAQP(obj)
            % Dimensions of the problem
            obj.m = size(obj.A,1);
            obj.n = size(obj.A,2);
            j = size(obj.P,1);

            if (j ~= obj.n)
                error('No valid matrix dimensions were introduced.');
            else
                obj.x = zeros(obj.n, obj.MaxIter+1);
                obj.z = zeros(obj.n, obj.MaxIter+1);
                obj.u = zeros(obj.m, 1);

                obj.rho = repmat(obj.rho, obj.n, 1);

                % Matrix equilibration
                [obj.Pt, obj.qt, obj.At, obj.c, obj.D, obj.E] = obj.Ruiz_equilibration(obj.P, obj.q, obj.A, obj.eps_equil);
            end
        end

        % General and eye solvers
        [x, z, Output] = quadratic_solve(obj);  
    end

    methods (Static, Access = private)
        [Pt, qt, At, c, D, E] = Ruiz_equilibration(P, q, A, eps);   % Matrix equilibration
    end
end