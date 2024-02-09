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
                % Feasible impulses, non-vertex
                U = [u -u];
                qf = [qf; qf];

                if (size(x,1) > 2 * n)
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, []); 
                        [x(2*n+1:end,:) ones(size(x,1)-2*n, 2 * n * N - size(x,2))]];
                else
                    X = [reshape(x(1:n,:), 1, []) reshape(x(n+1:2 * n,:), 1, [])];
                end
                Indx = X(1,:) ~= 0;
    
                if (xmax ~= Inf)
                    Indx = Indx & repmat( kron(x(2*n+3,:) < xmax, ones(1,n)), 1, 2);
                end
    
                if (xmin > 0)
                    Indx = Indx & repmat( kron(x(2*n+3,:) > xmin, ones(1,n)), 1, 2); 
                    Indx = Indx & repmat( kron(x(end,:) ~= 0, ones(1,n)), 1, 2);
                end
    
                U = U(:,Indx);                          % Considered subset of the sequence
                qf = qf(Indx,1).';                      % Cost function
                N = sum(Indx);
                
                % Assemble the constraint matrix 
                dim = size(x,1);
                if any(xmax ~= Inf)
                    U = [U zeros(m, 2*N); 
                        +eye(n*N) -eye(n*N) -kron([eye(N) zeros(N)],ones(n,1));
                        -eye(n*N) +eye(n*N) -kron([eye(N) zeros(N)],ones(n,1));
                         zeros(N, 2 * n * N) eye(N) eye(N)];                      % Complete matrix

                    lambda = [lambda; zeros(2 * n * N,1); xmax];            % Complete independent term
                    qf = [qf zeros(1,N)];                                   % Complete cost function
                end
    
                if any(xmin > 0)
                    if (dim > 2 * n)
                        U = [U zeros(size(U,1), N); eye(N) zeros(N) -eye(N)];
                    else
                        U = [U zeros(m, N); 
                            +eye(n*N) -eye(n*N) -kron(eye(N),ones(n,1));
                            -eye(n*N) +eye(n*N) -kron(eye(N),ones(n,1));
                            zeros(N, 2 * n * N) eye(N) -eye(N)];            % Complete matrix
                    end

                    lambda = [lambda; xmin];
                    qf = [qf zeros(1,N)];
                end

            end
            
        otherwise
    end

    if (null_flag)
       % Update of the SVD
        persistent svd_set
        persistent S
        persistent R
        persistent L

%         if isempty(svd_set)
%             [L, S, R] = svd(U);
%             svd_set = false;
%         else
%             % Rank 1 update 
%             R = [R; zeros(1,size(R,2))];
%             b = [zeros(1,size(R,2)) 1].';
%             [L1, S1, R1] = Solvers.Brand_SVD(L, S, R, U(:,end), b);
% 
%             [L, S, R] = svd(U);
%             norm(R-R1)
%             norm(L1-L)
%             norm(S1-S)
%         end

        % Compute the coordinate vector associated to a null impulse
        [L, S, R] = svd(U);
%         [Lb, Ub] = updateLU(U,2);
        e = U \ lambda;
        v = R(:, sum(S,1) == 0);
        v = v ./ sqrt(dot(v,v,1));
        
        c = ones(size(v,2),1);
        alpha = (e + v * c).';
        
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

%             if (any(mu == 0))
%                 idc = mu == 0;
%                 for i = 1:length(mu)
%                     % Rank 1 update 
%                     if (idc(i))
%                         b =  [repmat(0, 1, i-1) 1 repmat(0, 1, size(R,2)-i)];
%                         [L, S, R] = Solvers.Brand_SVD(L, S, R, -U(:,i), b.');
%                         L = L(:,~b);
%                         S = S(~b, ~b);
%                         R = R(~b,~b);
%                     end
%                 end
%             end
            
            % New sequence and slack variables
            switch (q)
                case 'L2'
                    V = mu .* V;
                    x(:,Indx) = reshape(V, size(x(:,Indx),2), size(x,1)).';
                                       
                case 'Linfty'
                    V = mu .* V;
                    X(:,Indx) = reshape(V, size(X(:,Indx),2), size(X,1)).';

                    if (size(x,1) > 2 * n)
                    else
                        dim = size(X,2) / 2;
                        x = [reshape(X(:,1:dim), n, []); reshape(X(:,dim+1:end), n, [])];
                    end
            end
        else
            warning('Infeasibility detected. Adding coasting will not improve the solution')
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
    y = L\b;
    x = U\y;
end