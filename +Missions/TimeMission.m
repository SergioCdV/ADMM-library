%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 09/10/25
% File: TimeMission.m 
% Issue: 0 
% Validated: 

%% Time optimal mission %% 
% Implementation of a class to define (possibly time-dependent) linear
% minimum time, fuel-optimal missions 

classdef TimeMission
  
    properties
        t;              % Mission clock/timeline discretization in the form of a row vector of dimensions 1 x N
        t0 = 0;         % Initial clock
        N = Inf;        % Maximum number of maneuvers in the mission
        Phi;            % State transition matrix evaluated at the mission clock. Matrix of dimensions m x mN, for a state vector in R^m
        B;              % Control input matrix. Matrix of dimensions m x nN, for a control vector in R^n
        x0;             % Initial conditions, as a vector of m x 1
        xf;             % Final conditions, as a vector of m x 1
    end

    properties (Hidden)
        m;              % State dimension
        n;              % Control dimension
    end

    methods
        % Constructor of the class
        function [obj] = TimeMission(tau, t0, myPhi, myB, myX0, myXF, myN)
            % Basic assignment
            obj.t = tau;        % Polynomial time
            obj.t0 = t0;        % Initial clock
            obj.Phi = myPhi;    % STM handle
            obj.B = myB;        % Input control matrix
            obj.x0 = myX0;      % Initial conditions
            obj.xf = myXF;      % Final conditions

            if (exist('myN', 'var'))
                if (myN > 0)
                    obj.N = myN;
                end
            end
        end
    end
end