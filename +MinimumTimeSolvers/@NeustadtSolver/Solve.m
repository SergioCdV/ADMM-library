%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: Solve.m 
% Issue: 0 
% Validated: 

%% Neustadt Solver %% 
% Main solver function %

% Inputs:  - object obj, the Linear Rendezvous Problem object
%          - scalar epsilon, numerical tolerance
%          - scalar rho, the augmented Lagrangian penalty parameter (rho > 0)
%          - scalar alpha, the overfitting parameter (2 > alpha > 0)

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, tf, obj] = Solve(obj, epsilon, rho, alpha)

    % Compute an initial guess for the minimum time and the control
    % magnitude

    % Continuate on the q-bound to reach the minimum time
    GoOn    = True;                       % Boolean to control the continuation process
    iter    = 1;                          % Initial iteration
    maxIter = 10;                         % Maximum number of iterations
    relTol  = 1E-9;                       % Relative assert tolerance
    bound_target = obj.Actuator.umax;     % Target q-norm bound

    tf           = zeros(1, maxIter);     % Pre-allocation of the final mission time

    while ( GoOn && iter < maxIter )
        % Update the STM and the mission constants 


        % Call the inner solver 
        [t, u, e, obj] = Solve@MinimumNormSolvers.NeustadtSolver(obj, rho, alpha);

        % Compute the bounds 
        qNorm = obj.Actuator.q.ComputeVectorNorm( u );

        % Convergence analysis
        if ( abs(qNorm - bound_target) <= relTol * max( qNorm, bound_target ) )
            GoOn = False;
            
        else
            % Update the guess on the final tf

            % Update all initial guesses 

            % Update the number of iterations 
            iter = iter + 1; 
        end
    end

    % Final results 
    tf = tf(1:iter);
end