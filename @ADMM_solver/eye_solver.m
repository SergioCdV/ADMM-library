        function [x, z, Output] = eye_solver(obj)
            % Initialization 
            iter = 1; 
            GoOn = true;

            tic;
            % Solver
            while (GoOn && iter <= obj.MaxIter)
                % X update
                obj.x(:,iter+1) = feval(obj.X_update, obj.x(:,iter), obj.z(:,iter), obj.u);

                % Z update with relaxation
                xh = obj.alpha * obj.x(:,iter+1) - (1-obj.alpha) * obj.z(:,iter);
                obj.z(:,iter+1) = feval(obj.Z_update, xh, obj.z(:,iter), obj.u);

                % Compute the residuals and update the 
                obj.u = obj.u + xh - obj.z(:,iter+1);
                
                % Convergence analysis 
                Output.objval(iter) = feval(obj.objective, obj.x(:,iter+1), obj.z(:,iter+1));
                Output.r_norm(iter) = norm(obj.x(:,iter+1) - obj.z(:,iter+1));
                Output.s_norm(iter) = norm(-obj.rho * (obj.z(:,iter+1) - obj.z(:,iter)));
            
                Output.eps_pri(iter) =  sqrt(obj.n) * obj.AbsTol + obj.RelTol * max([norm(obj.x(:,iter+1)), norm(obj.z(:,iter+1))]);
                Output.eps_dual(iter) = sqrt(obj.m) * obj.AbsTol + obj.RelTol * norm(obj.rho * obj.u);

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