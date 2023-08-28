%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 28/08/23
% File: thruster.m 
% Issue: 0 
% Validated: 

%% Thruster modelling %% 
% Implementation of a simple thruster model with its basic properties

classdef thruster
  
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
        function [obj] = thruster(myp, myumin, myumax)
            % Basic assignment
            obj.p = myp;

            switch (myp)
                case 'L1'
                    obj.q = 'Linfty';
                case 'L2'
                    obj.q = 'L2';
                case 'Linfty'
                    obj.q = 'L1';
                otherwise
                    error('No valid thruster configuration was selected');
            end

            if (exist('myumin', 'var'))
                obj.umin = myumin;
            end

            if (exist('myumax', 'var'))
                obj.umax = myumax;
            end
        end
    end
end