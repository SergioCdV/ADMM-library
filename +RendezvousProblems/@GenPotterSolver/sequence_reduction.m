function [x, cost, null_flag] = sequence_reduction(m, n, p, q, u, qf, x, lambda, xmax, xmin)
    % Final indices
    switch (q)
        case 'L2'
            Indx = x(1,:) ~= 0;                     % Feasible impulses, non-vertex      
            con = 0;
            
            % Add the constraints
            if (any(xmax ~= Inf))
                Indx = Indx & x(1,:) < xmax;
                Indx = Indx & x(2,:) > 0;
                con = con + 1;
            end

            if (any(xmin > 0))
                Indx = Indx & x(1,:) > xmin;
                Indx = Indx & x(end,:) > 0;
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
            con = 0;
    
            if (size(x,1) > 2 * n)
                X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2*n,:), 1, []); 
                    [reshape(x(2*n+1:3*n,:), 1, []) reshape(x(3*n+1:4*n,:), 1, [])]
                    ];

                if (size(x,1) > 4 * n)
                    X = [X; 
                        [reshape(x(4*n+1:5*n,:), 1, []) reshape(x(5*n+1:6*n,:), 1, [])]
                        ];
                end
            else
                X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, [])];
            end
    
            % NZ impulses
            Indx = X(1,:) ~= 0;
            len = size(Indx,2) / 2;
        
            if (any(xmax ~= Inf))
                ubnd_idx = X(2,:) > 0;
                ubnd_idx = ubnd_idx & [(+X(1,1:len) - X(1,len+1:end)) < repmat(xmax, 1, n) (-X(1,1:len) + X(1,len+1:end)) > -repmat(xmax, 1, n)];
                con = con + 1;
            end
    
            if (any(xmin > 0))
                if (any(xmax ~= Inf))
                    lbnd_idx = X(3,:) > 0;
                    lbnd_idx = lbnd_idx & [(+X(1,1:len) - X(1,len+1:end)) > repmat(xmin, 1, n) (-X(1,1:len) + X(1,len+1:end)) < -repmat(xmin, 1, n)];
                else
                    ubnd_idx = X(2,:) > 0;
                    ubnd_idx = ubnd_idx & [(+X(1,1:len) - X(1,len+1:end)) > repmat(xmin, 1, n) (-X(1,1:len) + X(1,len+1:end)) < -repmat(xmin, 1, n)];
                end
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
                    num_imp = (size(lambda,1) - m) / (2 + 2*n);
                                                  
                    % Pruning of constraints
                    con_idx = logical( kron(inf_idx, ones(1,n)) );

                    % Selection of variables
                    conv_idx = repmat(con_idx,1,2) & X(3,:) > 0;
                    idx1 = logical( [Indx inf_idx conv_idx inf_idx inf_idx] );
                    V = [X(1,Indx) X(2,inf_idx) X(3,conv_idx) X(4,inf_idx) X(5,inf_idx)];
                        
                    % Selection of equations
                    idx2 = logical( [ones(m,1); ...
                                     repmat([con_idx zeros(1,n*num_imp-size(X,2)/2)].', 2, 1); ...
                                     [inf_idx.'; zeros(num_imp-length(inf_idx),1)]; ...
                                     [inf_idx.'; zeros(num_imp-length(inf_idx),1)] ...
                                     ] );

                    U = u(idx2,idx1);                    % Considered subset of the sequence
                    qf = qf(idx1,1).';                   % Cost function
                    lambda = lambda(idx2,1);             % Independent term
    
                otherwise
                    % General selection of equations
                    num_imp = (size(lambda,1) - m) / (4 * n);

                    % Selection of variables
                    ubnd_idx = repmat(ubnd_idx(1:len) & ubnd_idx(len+1:end), 1, 2);
                    Indx = Indx & ubnd_idx;
                    idx1 = logical( [Indx ubnd_idx] );

                    V = [X(1,Indx) X(2,ubnd_idx)];
                       
                    % Selection of equations
                    idx2 = logical( [ones(m,1); ...
                                    repmat([ubnd_idx(1:len).'; zeros(n*num_imp-len,1)], 2, 1) ...
                                    ]);

                    U = u(idx2,idx1);                    % Considered subset of the sequence
                    qf = qf(idx1,1).';                   % Cost function
                    lambda = lambda(idx2,1);             % Independent term

                    cond_idx = Indx(1:len);
            end

            % Flag to indicate if the sequence is reducible
            cond_idx = reshape(Indx(1:len) | Indx(len+1:end), n, []);
            null_flag = sum( sum(cond_idx,1) > 0 ) > m;                
    end

    if (null_flag)
        % Update of the LU
        [~, Uf, ~] = lu(U);

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
                        len = size(X,2)/2;
                        dim = sum(ubnd_idx);

                        X(1,Indx) = V(1,1:sum(Indx));
                        X(2,ubnd_idx) = V(1,sum(Indx)+1:sum(Indx)+dim);
                        
                        x(1:n,:) = reshape(X(1,1:len), n, []);              % Positive side
                        x(n+1:2*n,:) = reshape(X(1,len+1:end), n, []);      % Negative side
                        x(2*n+1:3*n,:) = reshape(X(2,1:len), n, []);        % Infinity norm slacks
                        x(3*n+1:4*n,:) = reshape(X(2,len+1:end), n, []);    % Infinity norm slacks

                        if (size(x,1) == 6*n)
                            X(3,lbnd_idx) = V(1,sum(Indx)+dim+1:sum(Indx)+dim+sum(lbnd_idx));

                            x(4*n+1:5*n,:) = reshape(X(3,1:len), n, []);        % Infinity norm slacks
                            x(5*n+1:6*n,:) = reshape(X(3,len+1:end), n, []);    % Infinity norm slacks
                        end
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