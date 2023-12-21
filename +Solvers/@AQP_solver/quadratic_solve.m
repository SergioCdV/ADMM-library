%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 17/12/23
% File: quadratic_solve.m 
% Issue: 0 
% Validated: 

%% Quadratic Solve %%
% This function contains the implementation of an ADMM solver for general
% convex QP, based on the open-source OSQP

function [x, z, Output] = quadratic_solve(obj)
    % Initialization 
    iter = 1; 
    GoOn = true;            
    tic;

    % Solver
    while (GoOn && iter <= obj.MaxIter)
        % Compute the rho matrix 
        Rho = diag(obj.rho);
        index = logical( feval(obj.CheckCone, obj.z(:,iter)) );
        Rho(index,index) = Rho(index,index) * 1e3;

        % X update
        [xh, zh] = X_update(obj, Rho, obj.x(:,iter), obj.z(:,iter), obj.u);
        obj.x(:,iter+1) = obj.alpha * xh + (1 - obj.alpha) * obj.x(:,iter);

        % Z update with relaxation
        zh = obj.alpha * zh + (1-obj.alpha) * obj.z(:,iter) + obj.u ./ diag(Rho);
        obj.z(:,iter+1) = feval(obj.ConeProj, zh ./ diag(obj.E));

        % Compute the residuals
        obj.u = Rho * (zh - obj.z(:,iter+1));

        % Convergence analysis 
        Output.objval(iter) = objective(obj, obj.x(:,iter+1));
        Output.r_norm(iter) = norm(obj.At * obj.x(:,iter+1) - obj.z(:,iter+1), 'inf');
        Output.s_norm(iter) = norm(obj.Pt * obj.x(:,iter+1) + obj.At.' * obj.u + obj.qt, 'inf');
    
        a = max([norm(obj.At * obj.x(:,iter+1), 'inf'), norm(obj.z(:,iter+1), 'inf')]);
        Output.eps_pri(iter) = obj.AbsTol + obj.RelTol * a;

        b = max([norm(obj.Pt * obj.x(:,iter+1), "inf"), norm(obj.At.' * obj.u, "inf"), norm(obj.qt,'inf')]);
        Output.eps_dual(iter) = obj.AbsTol + obj.RelTol * b;

        % Update the penalty matrix 
        if (a ~= 0)
            a = Output.r_norm(iter) / a;              % Primal ratio
            if (b ~= 0)
                b = Output.eps_dual(iter) / b;        % Dual ratio
                obj.rho = obj.rho * sqrt(a / b);      % Penalty factors
            end
        end
        
        if (obj.QUIET)
            if (iter == 1)
                fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', 'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
                fprintf('-----------------------------------------------------------------------------------\n');
                fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter), Output.objval(iter));
            else
                fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, Output.r_norm(iter), Output.eps_pri(iter), Output.s_norm(iter), Output.eps_dual(iter), Output.objval(iter));
            end
        end

        % Update the state machine
        if (Output.r_norm(iter) < Output.eps_pri(iter) && Output.s_norm(iter) < Output.eps_dual(iter))
            GoOn = false;
        else
            iter = iter + 1;
        end
    end

    % Final results 
    Output.Time = toc;
    Output.Result = ~GoOn; 
    
    if (GoOn)
        x = obj.x(:,1:iter-1); 
        z = obj.z(:,1:iter-1);
        Output.Iterations = iter-1;
    else
        x = obj.x(:,1:iter); 
        z = obj.z(:,1:iter);
        Output.Iterations = iter;
    end
end

%% Auxiliary functions 
% Update on the splitting X variables (QP problem) 
function [xt, zt] = X_update(obj, Rho, x, z, y)
    % KKT matrix 
    A = obj.Pt + obj.sigma * eye(size(obj.Pt)) + (obj.At.' * Rho * obj.At);
    b = obj.sigma * x - obj.qt + obj.At.' * (Rho * z - y);

    % Final variables
    xt = A \ b; 
    zt = obj.At * xt;
end

% Objective function 
function [obj] = objective(obj, x)
    obj = 0.5 * x.' * obj.Pt * x + obj.qt.' * x; 
end
