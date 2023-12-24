%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: NeustadtSolver.m 
% Issue: 0 
% Validated: 

%% Neustadt Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

classdef QuadNeustadtSolver < RendezvousProblems.SolverRendezvous
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
        function [obj] = QuadNeustadtSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj, rho, alpha);
    end

    % ADMM functions
    methods (Static)
        [x] = ConeCheck(indices, q, b, z);
        [z] = ConeProj(indices, q, b, z);
    end

end