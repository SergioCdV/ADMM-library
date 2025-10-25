%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: ProblemTypes.m 
% Issue: 0 
% Validated: 

%% Problem Types %%
% Enumerate for different types of boundary value problems

classdef ProblemTypes < int8
   enumeration
       IVP     (1)
       TBVP    (2)
   end
end
