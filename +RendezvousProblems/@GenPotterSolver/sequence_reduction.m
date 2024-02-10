function [x, cost, null_flag] = sequence_reduction(m, n, p, q, u, qf, x, lambda, xmax, xmin)
    % Final indices
    switch (q)
        case 'L2'
            Indx = x(1,:) ~= 0;                     % Feasible impulses, non-vertex      
            con = 0;
            
            % Add the constraints
            if (any(xmax ~= Inf))
                Indx = Indx & x(1,:) < xmax;
                Indx = Indx & x(2,:) ~= 0;
                con = con + 1;
            end

            if (any(xmin > 0))
                Indx = Indx & x(1,:) > xmin;
                Indx = Indx & x(end,:) ~= 0;
                con = con + 2;
            end

            N = sum(Indx);                          % Number of non-zero variables
            V = reshape(x(:,Indx).', 1, []);        % Get all the variables in a row
            null_flag = N > m;                      % Flag to indicate if the sequence is reducible

            % Reduce the linear problem
            switch (con)
                case 0
                    % General selection of equations
                    U = u(:,Indx);                       % Considered subset of the sequence
                    qf = qf(Indx,1).';                   % Cost function

                case 3
                    % General selection of equations
                    num_imp = ( size(lambda,1) - m ) / 2;
                    idx1 = logical( repmat(Indx, 1, 3) );    
                    idx2 = logical( [ones(m,1); repmat( [Indx.'; zeros(num_imp-size(x,2),1)], 2, 1) ] );

                    U = u(idx2,idx1);                    % Considered subset of the sequence
                    qf = qf(idx1,1).';                   % Cost function
                    lambda = lambda(idx2,1);             % Independent term

                otherwise
                    % General selection of equations
                    num_imp = size(lambda,1) - m;
                    idx1 = logical( repmat(Indx, 1, 2) );    
                    idx2 = logical( [ones(m,1); Indx.'; zeros(num_imp-size(x,2),1)] );

                    U = u(idx2,idx1);                    % Considered subset of the sequence
                    qf = qf(idx1,1).';                   % Cost function
                    lambda = lambda(idx2,1);             % Independent term
            end

        case 'Linfty'
            Indx = sum(x(1:2*n,:)) ~= 0;
            N = sum(Indx);                          % Number of non-zero impulses
            null_flag = N > m;                      % Flag to indicate if the sequence is reducible

            if (null_flag)    
                con = 0;
                N = size(x,2);

                if (size(x,1) > 2 * n)
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, []); 
                        [x(2*n+1:end,:) ones(size(x,1)-2*n, (2 * n -1) * N)]
                        ];
                else
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, [])];
                end
                Indx = X(1,:) ~= 0;
    
                if (any(xmax ~= Inf))
                    Indx = Indx & repmat( kron(x(2*n+1,:) < xmax, ones(1,n)), 1, 2);
                    Indx = Indx & repmat( kron(x(2*n+2,:) ~= 0, ones(1,n)), 1, 2);
                    con = con + 1;
                end
    
                if (any(xmin > 0))
                    Indx = Indx & repmat( kron(x(2*n+1,:) > xmin, ones(1,n)), 1, 2); 
                    Indx = Indx & repmat( kron(x(end,:) ~= 0, ones(1,n)), 1, 2);
                    con = con + 2;
                end

                % Reduce the linear problem
                switch (con)
                    case 0
                        % General selection of equations
                        U = u(:,Indx);                       % Considered subset of the sequence
                        qf = qf(Indx,1).';                   % Cost function
                        V = X(1,Indx);                       % Get all the variables in a row
    
                    case 3
                        % General selection of equations
                        num_imp = ( size(lambda,1) - m ) / 2;
                        idx1 = logical( repmat(Indx, 1, 3) );    
                        idx2 = logical( [ones(m,1); repmat( [Indx.'; zeros(num_imp-N,1)], 2, 1) ] );
    
                        U = u(idx2,idx1);                    % Considered subset of the sequence
                        qf = qf(idx1,1).';                   % Cost function
                        lambda = lambda(idx2,1);             % Independent term
    
                    otherwise
                        % General selection of equations
                        num_imp = (size(lambda,1) - m) / (1 + 2*n);

                        inf_idx = logical( sum( [reshape(X(1,1:n*N), n, N); reshape(-X(1,n*N+1:end), n, N)], 1) );
                        idx1 = logical( [Indx repmat(inf_idx, 1, 2)] );   

                        idx2 = logical( [ones(m,1); ...
                                         [Indx.'; zeros(2*n*num_imp-length(Indx),1)]; ...
                                         [inf_idx.'; zeros(num_imp-length(inf_idx),1)] ...
                                         ] );
    
                        U = u(idx2,idx1);                    % Considered subset of the sequence
                        qf = qf(idx1,1).';                   % Cost function
                        lambda = lambda(idx2,1);             % Independent term

                        V = [X(1,Indx) X(2,inf_idx) X(3,inf_idx)];
                end
            end
            
        otherwise
    end

    if (null_flag)
       % Update of the LU
       [~, ~, P] = lu(U, 'vector');
       [~, Uf, ~] = lu(U(P,:), 'vector');
%         [Lb, Ub] = updateLU(U,2);

        % Compute the coordinate vector associated to a null impulse
        M = size(Uf,1);
        d = size(U,2) - M;
        v = [zeros(M,d); eye(d)];
        if (det(Uf(:,1:M)) ~= 0)
            for i = 1:d
                v(1:M,i) = BTRAN(-Uf(:,1:M), Uf(:,M+i));        % Basis for the nullspace 
            end
        else
            v = null(Uf);
        end
        
        c = ones(size(v,2),1);
        alpha = (v * c).';
        
        if (sum(lambda) == 0)
            Alpha = dot( qf, alpha );           % New cost function
            if (Alpha < 0)
                alpha = -alpha;                 % Ensure feasibility of the solution in the unconstrained case
            end
        end
    
        % Update the impulse sequence 
        beta = V ./ alpha;                      % Compute the sequence ratios
        beta_r = min(beta(alpha > 0));          % Minimum ratio

        if (beta_r > 0)
            mu = 1 - beta_r ./ beta;            % New sequence magnitudes
            
            % New sequence and slack variables
            switch (q)
                case 'L2'
                    V = mu .* V;
                    x(:,Indx) = reshape(V, size(x(:,Indx),2), size(x,1)).';
                                       
                case 'Linfty'
                    V = mu .* V;
                    
                    if (size(x,1) > 2 * n)
                        dim = sum(inf_idx);
                        X(1,Indx) = V(1,1:sum(Indx));
                        X(2:end,inf_idx) = reshape(V(1,sum(Indx)+1:end), dim, size(X(2:end,:),1)).';

                        len = size(X,2)/2;
                        x = [reshape(X(1,1:len), n, []);               % Positive side
                             reshape(X(1,len+1:end), n, []);           % Negative side
                             X(2:end,1:length(inf_idx))                % Infinity norm and slacks
                             ];
                    else
                        dim = size(X,2) / 2;
                        X(:,Indx) = reshape(V, size(X(:,Indx),2), size(X,1)).';
                        x = [reshape(X(:,1:dim), n, []); reshape(X(:,dim+1:end), n, [])];
                    end
            end
        else
            warning('Infeasibility detected. Adding coasting will not improve the solution')
            null_flag = false;                      % Flag to indicate if the sequence is reducible
        end
    end

    % Final cost and saturation
    switch (p)
        case 'L2'
            cost = x(1,:);
            null_flag = false;                      % Flag to indicate if the sequence is reducible
            
        case 'L1'
            cost = sum( abs(x(1:2*n,:)), 1 );

        case 'Linfty'
            cost = max( abs(dV(1:n,:)), [], 1 );
    end

    % Final cost
    cost = sum(cost);
end

%% Auxiliary functions 
% Solve a linear system using LU
function [x] = LUsolve(L,U,b)
    y = FTRAN(L,b);
    x = BTRAN(U,y);
end

% Forward substitution
function [x] = FTRAN(L, b)

    x = zeros(size(L,2),1);
    for i = 1:size(L,2)
        x(i,1) = ( b(i,1) - dot(L(i,1:i-1), x(1:i-1,1)) ) / L(i,i);
    end

end

% Backward substitution
function [x] = BTRAN(U, b)
    
    n = size(U,2);
    x = zeros(size(U,2),1);
    for i = n:-1:1
        x(i,1) = ( b(i,1) - dot(U(i,i+1:n), x(i+1:n,1)) ) / U(i,i);
    end

end