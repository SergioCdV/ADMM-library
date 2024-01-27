%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 27/01/24
% File: GenPotterSlver.m 
% Issue: 0 
% Validated: 

%% Generalized Potter Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

classdef GenPotterSolver < RendezvousProblems.SolverRendezvous
    % Basic properties
    properties
        t;                          % Execution clocks
        u;                          % Impulsive control law
        e;                          % Final missvector
        Cost;                       % Control law cost
        SolveTime;                  % Elapsed time 

        pol_flag = false;           % Boolean flag to activate polishing
    end

    % Methods
    methods
        % Constructor
        function [obj] = GenPotterSolver(myMission, myThruster)
            obj@RendezvousProblems.SolverRendezvous(myMission, myThruster);
        end

        % Solver
        [t, u, e, obj] = Solve(obj);
    end

    methods (Static)
        [dV, cost] = PVT_pruner(Phi, B, dV, dVmin, dVmax, p);
    end

end