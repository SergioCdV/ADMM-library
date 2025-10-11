%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: FuelMission.m 
% Issue: 0 
% Validated: 

%% Fuel optimal mission %% 
% Implementation of a class to define (possibly time-dependent) linear
% minimum norm, fuel-optimal missions 

classdef FuelMission
  
    properties
        t;              % Mission clock/timeline discretization in the form of a row vector of dimensions 1 x N
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
        function [obj] = FuelMission(t, myPhi, myB, myX0, myXF, myN)
            % Basic assignment
            obj.t = t;
            obj.Phi = myPhi;
            obj.B = myB;
            obj.x0 = myX0; 
            obj.xf = myXF;

            if (exist('myN', 'var'))
                if (myN > 0)
                    obj.N = myN;
                end
            end

            % Timeline check
            dt = diff(t); 
            if (any(dt < 0))
                error('The mission timeline is not strictly positive');
            end

            % Dimensionality check 
            if (size(obj.B,1) ~= size(obj.Phi,1))
                error('State space model is correctly defined');
            end

            if (length(t) * size(obj.Phi,1) ~= size(obj.Phi,2))
                error('The STM is not completely provided over all the mission clock');
            else
                obj.m = size(obj.Phi,1);
            end

            if (mod(size(obj.B,2),length(t)) ~= 0)
                error('The control input matrix is not completely provided over all the mission clock');
            else
                obj.n = size(obj.B,2) / length(t);
            end
        end
    end
end