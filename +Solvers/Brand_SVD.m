%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 06/02/24
% File: Brand_SVD.m 
% Issue: 0 
% Validated: 

%% Incremental SVD %%
% This function contains the implementation of Brand's method for rank-1
% updates of SVD decompositions

function [U, S, V] = Brand_SVD(U, S, V, a, b)
    % Initilization 
    r = size(U,2);
    
    m = U.' * a; 
    c = U * m;
    p = a - c; 
    
    p_norm = sqrt( dot(p,p) );
    if (p_norm > 0)
        P = orth(p);
    else
        P = p;
    end

    n = V.' * b; 
    d = V * n;
    q = b - d; 
    q_norm = sqrt( dot(q,q) );

    if (q_norm > 0)
        Q = orth(q);
    else
        Q = q;
    end

    % Main process
    K = [S zeros(size(S,1),1); zeros(1,size(m,1)) 0] + [m * n.' q_norm * m; p_norm * n.' q_norm * p_norm];
    [Up, S, Vp] = svds(K,r);
    
    U = [U P] * Up;
    V = [V Q] * Vp;

    if (p_norm == 0)
        U = U(:,1:end-1);
    end

    if (q_norm == 0)
        V = V(:,1:end-1);
    end
end