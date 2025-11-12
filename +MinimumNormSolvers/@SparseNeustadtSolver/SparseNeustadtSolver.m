%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 11/11/25
% File: SparseNeustadtSolver.m 
% Issue: 0 
% Validated: 

%% Sparse Neustadt Solver %% 
% Implementation of a object-oriented solver for minimum norm problems via ADMM

classdef SparseNeustadtSolver < MinimumNormSolvers.SolverMinimumNorm
    % Basic properties
    properties
        t;                  % Execution clocks
        u;                  % Impulsive control law
        e;                  % Final missvector
        Cost;               % Control law cost
        SolveTime;          % Elapsed time 
        Report;             % Optimization report
    end

    % Methods
    methods
        % Constructor
        function [obj] = SparseNeustadtSolver(myMission, myActuator)
            obj@MinimumNormSolvers.SolverMinimumNorm(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, epsilon, rho, alpha, init_guess);
    end

    methods (Static)
        % ADMM functions
        [x] = x_update(F, V, c, rho, x, z, u);
        [z] = z_update(m, n, nx, q, c, b, F, rho, x, z, u);
    end

end