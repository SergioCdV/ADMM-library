%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: PrimalSolver.m 
% Issue: 0 
% Validated: 

%% Primal Solver %% 
% Implementation of a object-oriented solver for minimum norm problems via ADMM

classdef PrimalSolver < MinimumNormSolvers.SolverMinimumNorm
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
        function [obj] = PrimalSolver(myMission, myActuator)
            obj@MinimumNormSolvers.SolverMinimumNorm(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha, equil_flag);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(p, x);
        [x] = x_update(pInvA, Atb, x, z, u);
        [z] = z_update(indices, p, q, umin, umax, K, rho, x, z, u);
    end
end