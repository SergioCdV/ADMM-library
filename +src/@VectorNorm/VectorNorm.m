%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: VectorNorm.m 
% Issue: 0 
% Validated: 

%% Vector norms %%
% This class definition provides an enumeration for the different norms to
% characterize fuel consumption

classdef VectorNorm < double
   enumeration
       L1     (1)
       L2     (2)
       Linfty (Inf)
   end

   methods 
       % Compute the norm associated to a given enumeration value and input vector
       [norm] = ComputeVectorNorm( obj, myVector );

       % Compute a Holders' conjugate
       [conj] = HolderConjugate(obj);
   end

   methods (Static)
       % Get the enumeration associated to a given double
       [enum] = GetEnum( myNum );
   end
end