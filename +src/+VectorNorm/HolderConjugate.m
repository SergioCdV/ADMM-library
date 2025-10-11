%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: HolderConjugate.m 
% Issue: 0 
% Validated: 

%% Holder's Conjugate %%
% This function computes the Holder conjugate to a given vector norm

function [conj] = HolderConjugate( obj )
    conj = 1 / ( 1 - double(obj) ); 
    conj = obj.GetEnum( conj );
end