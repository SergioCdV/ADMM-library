%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: GetEnum.m 
% Issue: 0 
% Validated: 

%% Get enumeration %%
% This function transforms a given double into a valid VectorNorm object

function [enum] = GetEnum( myNum )
    
    allEnums = enumeration('VectorNorm');     
    allVals  = double(allEnums);         

    idx = find(allVals == myNum, 1, 'first');
    
    if isempty(idx)
        error('The input argument cannot be converted into a supported vector norm');
    else
        enum = allEnums(idx);
    end
end