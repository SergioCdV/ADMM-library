%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Primal_solver.m 
% Issue: 0 
% Validated: 

%% Rendezvous Solver %% 
% Implementation of a object-oriented, abstract solver for rendezvous problems
% via ADMM

classdef SolverRendezvous
    % Basic properties
    properties
        Mission;            % Mission to be solved 
        Thruster;           % Thruster to be used
    end

    % Fundamental methods
    methods
        % Constructor
        function [obj] = SolverRendezvous(myMission, myThruster)
            obj.Mission = myMission;
            obj.Thruster = myThruster;
        end

        % Solver
        [t, u, e, obj] = Solve(obj);
    end

end