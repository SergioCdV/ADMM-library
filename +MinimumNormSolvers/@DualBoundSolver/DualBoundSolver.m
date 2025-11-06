%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: DualBoundSolver.m 
% Issue: 0 
% Validated: 

%% Neustadt Solver %% 
% Implementation of a object-oriented solver for minimum norm problems via ADMM

classdef DualBoundSolver < MinimumNormSolvers.NeustadtSolver
    properties (Access = private)
        window_ratio = 1.0;            % Percentage of time windows over the grid
    end

    % Methods
    methods
        % Constructor
        function [obj] = DualBoundSolver(myMission, myActuator)
            obj@MinimumNormSolvers.NeustadtSolver(myMission, myActuator);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, epsilon, rho, alpha, init_guess);
    end

    % ADMM functions
    methods (Static)
        [x] = x_update(Phi, c, idx, rho, x, z, u);
        [z] = z_update(n, sigma_pos, sigma_map, sigma_unique, q, Phi, b, ~, x, ~, u);
    end

end