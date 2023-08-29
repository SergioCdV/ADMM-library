%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 29/08/23
% File: PrimalSolver.m 
% Issue: 0 
% Validated: 

%% Carter Rendezvous Solver %% 
% Implementation of a object-oriented solver for dual linear rendezvous problems
% via ADMM

classdef CarterSolver < RendezvousProblems.SolverRendezvous
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
        function [obj] = CarterSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [p] = objective(p, cum_part, z);
        [x] = x_update(n, q, pInvA, Atb, x, z, u);
        [z] = z_update(indices, p, q, umin, umax, K, Phi, b, rho, x, z, u);
    end

end