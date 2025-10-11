%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: Actuator.m 
% Issue: 0 
% Validated: 

%% Actuator modelling %% 
% Implementation of a simple actuator model with its basic properties

classdef Actuator
  
    properties
        p;              % Fuel consumption norm
        umin = 0;       % Minimum control authority
        umax = Inf;     % Maximum control authority
    end

    properties (Hidden)
        q;              % Dual control authority norm
    end

    methods
        % Constructor of the class
        function [obj] = Actuator(myp, myumin, myumax)
            % Basic assignment
            obj.p = myp;
            obj.q = obj.p.HolderConjugate();

            if (exist('myumin', 'var'))
                obj.umin = myumin;
            end

            if (exist('myumax', 'var'))
                obj.umax = myumax;
            end
        end
    end
end