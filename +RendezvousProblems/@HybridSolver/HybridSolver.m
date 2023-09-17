%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: PrimalSolver.m 
% Issue: 0 
% Validated: 

%% Hybrid Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

classdef HybridSolver < RendezvousProblems.SolverRendezvous
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
        function [obj] = HybridSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(cum_part, x, z);
        [x] = x_update(Phi, b, rho, x, z, u);
        [z] = z_update(indices, umin, umax, K, rho, x, z, u);
    end

end