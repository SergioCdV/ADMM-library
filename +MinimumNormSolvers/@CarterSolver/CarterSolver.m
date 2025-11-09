%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 29/08/23
% File: CarterSolver.m 
% Issue: 0 
% Validated: 

%% Carter Solver %% 
% Implementation of a object-oriented solver for minimum norm problems via ADMM

classdef CarterSolver < MinimumNormSolvers.SolverMinimumNorm
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
        function [obj] = CarterSolver(myMission, myActuator)
            obj@MinimumNormSolvers.SolverMinimumNorm(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha, equil_flag);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(p, z);
        [x] = x_update(n, q, pInvA, Atb, x, z, u);
        [z] = z_update(n, p, q, umin, umax, K, Phi, b, rho, x, z, u);
    end

end