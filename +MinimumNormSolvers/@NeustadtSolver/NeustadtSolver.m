%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: NeustadtSolver.m 
% Issue: 0 
% Validated: 

%% Neustadt Solver %% 
% Implementation of a object-oriented solver for minimum norm problems via ADMM

classdef NeustadtSolver < MinimumNormSolvers.SolverMinimumNorm
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
            obj@MinimumNormSolvers.SolverMinimumNorm(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, epsilon, rho, alpha, init_guess);
    end

    methods (Static)
        % ADMM functions
        [p] = objective(c, x);
        [x] = x_update(Phi, c, rho, x, z, u);
        [z] = z_update(n, q, Phi, b, rho, x, z, u);
        
        % Input reconstruction
        [t, u] = ImpulseReconstruction(t, b, Y, p_norm, a, epsilon)
    end

end