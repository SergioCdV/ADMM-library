%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: PrimalSolver.m 
% Issue: 0 
% Validated: 

%% Primal Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

classdef GPrimalSolver < RendezvousProblems.SolverRendezvous
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
        function [obj] = GPrimalSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(p, cum_part, x);
        [x] = x_update(indices, p, rho, x, z, u);
        [z] = z_update(indices, p, q, umin, umax, K, pInvA, Atb, rho, x, z, u);
    end

end