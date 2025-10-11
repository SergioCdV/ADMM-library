%% Optimal control via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: VectorNorm.m 
% Issue: 0 
% Validated: 

%% Compute Vector Norm %%
% This functions takes the norm of myVector as specificed by myNorm in
% a vectorized fashion

function [norm] = ComputeVectorNorm( obj, myVector )
    % Compute the norm
    switch (obj)
        case src.VectorNorm.L2
            norm = sqrt( dot(myVector, myVector, 1) );
            
        case src.VectorNorm.L1
            norm = sum( abs(myVector), 1 );

        case src.VectorNorm.Linfty
            norm = max( abs(myVector), [], 1 );
    end
end