%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: System.m 
% Issue: 0 
% Validated: 

%% System %%
% This class presents an abstract implementation of dynamical systems

classdef System 
    properties
        Dynamics;          % First-order vector field of the evolution
        ProblemType;       % Type of IVP to be solved
        BC;                % Boundary conditions
    end

    methods
        % Basic constructor
        function [obj] = System(myField, myProblemType, myBC)
            obj.Dynamics    = myField; 
            obj.ProblemType = myProblemType;
            obj.BC = myBC;
        end
    end

end