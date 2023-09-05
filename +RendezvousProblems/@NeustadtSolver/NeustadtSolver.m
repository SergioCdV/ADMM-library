%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: NeustadtSolver.m 
% Issue: 0 
% Validated: 

%% Neustadt Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

classdef NeustadtSolver < RendezvousProblems.SolverRendezvous
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
        function [obj] = NeustadtSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(c, x);
        [x] = x_update(m, Phi, c, Phi2, rho, x, z, u);
        [z] = z_update(indices, q, K, Phi, b, rho, x, z, u);
    end

end