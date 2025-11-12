%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 09/11/25
% File: DecompositionUV.m 
% Issue: 0 
% Validated: 

%% U-V decomposition %%
% This function contains the implementation of the UV decomposition for
% sparse matrices in 'A First-Order Numerical Algorithm without Matrix
% Operations'

function [U, V] = DecompositionUV(A)
    % Compute the number of non-zero elements
    [m, n] = size(A); 

    pos = A ~= 0;
    N      = sum( pos, 1 );
    Ncum   = cumsum( [1 N] );
    sparse = sum( N, 2 );

    if ( sparse )
        % Perform the decomposition 
        U = zeros(m, sparse);
        V = zeros(n, sparse);
    
        % Compute the decomposition 
        for i = 1:n
            rows_i   = find( pos(:,i) );
            cols_i   = Ncum(i): Ncum(i+1)-1;
            values_i = A(rows_i, i);
            
            % V matrix
            V(i, cols_i) = 1;            

            % U matrix
            U(rows_i, cols_i) = diag( values_i );
        end
    else
        U = zeros(m, n); 
        V = zeros(n, n);
    end
end