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
    m = U.' * a; 
    c = U * m;
    p = a - c; 
    p_norm = norm(p);

    if (p_norm > 0)
        P = p / p_norm;
        p_flag = true;
    else
        P = zeros(size(p));
        p_flag = false;
    end

    n = V.' * b; 
    d = V * n;
    q = b - d; 
    q_norm = norm(q);

    if (q_norm > 0)
        Q = q / q_norm;
        q_flag = true;
    else
        Q = zeros(size(q));
        q_flag = false;
    end

    % Main process
    K = [S zeros(size(S,1),1); zeros(1,size(S,2)) 0] + [m; p_norm]*[n; q_norm].';
    [Up, S, Vp] = svd(K, 'econ');
    
    % Final update
    if (1)
        U = [U P] * Up; 
    else
        U = U * Up(1:end-1,1:end-1);
        S = S(1:1:end-1,:);
    end

    if (1)
        V = [V Q] * Vp; 
    else
        V = Vp(:,1:end-1);
        S = S(:,1:end-1);
    end
    
end