%% Modified Chebyshev Picard Iteration Solver %% 
% Sergio Cuevas del Valle
% Date: 19/10/25
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Solve %%
% This function implements the iteration procedure to solve a given BP

function [time, x, Output] = Solve(obj, Problem, x0, delta_t)
    % Pre-allocation 
    m = size(x0,1);                         % Dimension of the state space
    X = zeros(m * obj.N, obj.maxIter);      % Evolution of the trajectory
    X0 = Problem.BC;                        % Problem boundary conditions     

    % Set up of the method 
    GoOn = true;                            % Convergence boolean
    iter = 1;                               % Initial iteration
    error = 1E3;                            % Initial error
    
    switch Problem.ProblemType
        case Systems.ProblemTypes.IVP
            hndl_ = @(x, beta) ( beta(1:2,:) + [2 * x(1,:); zeros(1, size(x,2))] );

        case Systems.ProblemTypes.TBVP
            hndl_ = @(x, beta) ( [  x(1,:) + x(end,:) - 2 * sum( beta(3:2:end,:), 1 ); ...
                                    0.5 * ( x(end,:) - x(1,:) ) - sum( beta(4:2:end,:), 1 ) ...
                                 ] );
    end

    % Dynamical system handle 
    [time, J] = obj.PhysicalTime( delta_t );
    dynamics = @(s)( delta_t * J .* Problem.Dynamics(time, s) );

    % Main loop
    while ( GoOn && iter < obj.maxIter )
        % Initial iteration 
        X(:,iter) = reshape( x0, [], 1 );

        % Evaluate the dynamics 
        g = dynamics( x0 ).'; 

        % Compute the beta coefficients as a function of the problem
        beta = obj.Ca * g;                                              % Inner CGL nodes coefficients
        beta(1:2,:) = hndl_( X0.', beta );                              % Adapt first coefficients depending on the problem type

        % Compute the new state trajectory
        x = beta.' * obj.Cx.';

        % Convergence check 
        dX = x - x0; 
        dX = sqrt( dot(dX, dX, 2) );

        if ( error < obj.tol && norm(dX) < obj.tol )
            GoOn = false;
        else
            % Update the iteration scheme
            x0 = x;
            error = norm(dX);
            iter = iter + 1;
        end
    end

    % Final output 
    x = x.';

    Output.Result     = ~GoOn; 
    Output.Iterations = iter; 
    Output.Error      = error; 
    Output.Evolution  = X(:,1:iter);
end