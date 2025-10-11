%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: NeustadtSolver.m 
% Issue: 0 
% Validated: 

%% Neustadt Solver %% 
% Implementation of a object-oriented solver for minimum time problems via ADMM

classdef NeustadtSolver < MinimumNormSolvers.NeustadtSolver
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
        function [obj] = NeustadtSolver(myMission, myActuator)
            obj@MinimumNormSolvers.NeustadtSolver(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(c, x);
        [x] = x_update(m, Phi, c, rho, x, z, u);
        [z] = z_update(indices, q, Phi, b, rho, x, z, u);
    end
end