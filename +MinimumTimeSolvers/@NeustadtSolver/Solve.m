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

function [t, u, e, tf, obj] = Solve(obj, epsilon, rho, alpha, beta)
    % Sanity checks 
    if ( ~exist('beta', 'var') ) 
        beta = 0.01;
    end

    if (~exist('alpha', 'var'))
        alpha = 1;
    end

    % Continuate on the q-bound to reach the minimum time
    GoOn    = true;                       % Boolean to control the continuation process
    iter    = 1;                          % Initial iteration
    maxIter = 100;                         % Maximum number of iterations
    relTol  = 1E-9;                       % Relative assert tolerance
    bound_target = obj.Actuator.umax;     % Target q-norm bound
    
    tf      = zeros(1, maxIter);          % Pre-allocation of the final mission time

    % Compute an initial guess for the minimum time 
    tf(1) = obj.Mission.t0 + ( norm(obj.Mission.xf) - norm(obj.Mission.x0) ) / obj.Actuator.umax;

    while ( GoOn && iter < maxIter )
        % Update the STM and the mission constants 
        t   = ( tf(iter) - obj.Mission.t0 ) * obj.Mission.t;
        Phi = obj.Mission.Phi( tf(iter) );

        fuelMission = Missions.FuelMission( t, Phi, obj.Mission.B, obj.Mission.x0, obj.Mission.xf, obj.Mission.N );

        % Call the inner solver 
        NormSolver   = MinimumNormSolvers.NeustadtSolver(fuelMission, obj.Actuator);
        [t, u, e, ~] = NormSolver.Solve(epsilon, rho, alpha);

        % Compute the bounds 
        qNorm = obj.Actuator.q.ComputeVectorNorm( u );

        % Convergence analysis
        if ( abs(qNorm - bound_target) <= relTol * max( qNorm, bound_target ) )
            GoOn = False;
            
        else
            % Update the guess on the final tf
            ratio = qNorm / bound_target;
            tf(iter + 1) = tf(iter) + min(pi, max(obj.Mission.t0, tf(iter) + beta * ratio))
%             beta = beta * ratio;

            % Update all initial guesses 

            % Update the number of iterations 
            iter = iter + 1; 
        end
    end

    % Final results 
    tf = tf(1:iter);
end