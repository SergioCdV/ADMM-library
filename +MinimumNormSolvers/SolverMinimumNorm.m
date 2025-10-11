%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: SolverMinimumNorm.m 
% Issue: 0 
% Validated: 

%% Minimum Norm Solver %% 
% Implementation of a object-oriented, abstract solver for minimum norm problems via ADMM

classdef SolverMinimumNorm
    % Basic properties
    properties
        Mission;            % Mission to be solved 
        Actuator;           % Actuator to be used
    end

    % Fundamental methods
    methods
        % Constructor
        function [obj] = SolverMinimumNorm(myMission, myActuator)
            obj.Mission  = myMission;
            obj.Actuator = myActuator;
        end

        % Solver
        [t, u, e, obj] = Solve(obj);
    end

end