%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 26/11/25
% File: GivensRotation.m 
% Issue: 0 
% Validated:

%% Givens rotation %%
% This function performs Givens rotations over a given matrix (right multiplication), 
% as described in https://normalsplines.blogspot.com/2019/02/algorithms-for-updating-cholesky.html

function [L] = GivensRotation( A, r )
    % Pre-allocation
    B = A;                
    [n, m] = size(A);  

    for t = 0:m-r-1
        % Compute Givens parameter 
        k = r + t;
        k = min(k, n);
        b = B(k,k+1);            % l_{r+t, r+t+1}

        if ( b ~= 0 )
            a = B(k,k);           % l_{r+t, r+t}

            % Numericall-stable computation of the rotation
            denom = hypot(a, b);  
            c = a / denom;
            s = b / denom;

            % Propagate the rotation
            for i = k:n
                old1 = B(i,k);
                old2 = B(i,k+1);
                B(i,k)   =  c * old1 + s * old2;
                B(i,k+1) = -s * old1 + c * old2;
            end
        end
    end
    
    % Final output
    L = B;
end