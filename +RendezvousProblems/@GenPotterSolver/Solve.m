%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 27/01/24
% File: Solve.m 
% Issue: 0 
% Validated: 

%% General Potter Rendezvous Solver %% 
% Implementation of a object-oriented solver for primal linear rendezvous problems
% via ADMM

% Inputs:  - object obj, the Linear Rendezvous Problem object

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)
%          - vector e, of dimensions m x 1, the final rendezvous missvector

function [t, u, e, obj] = Solve(obj)
    % Preallocation 
    x0 = obj.Mission.x0;                    % Initial conditions 
    xf = obj.Mission.xf;                    % Final conditions 

    t = obj.Mission.t;                      % Mission clock
    N = length(t);                          % Number of total opportunities
    STM = obj.Mission.Phi;                  % STM of the system
    B = obj.Mission.B;                      % Control input of the system
    m = obj.Mission.m;                      % State vector dimension
    n = obj.Mission.n;                      % Control input dimension

    Phi = zeros(size(STM,1), size(B,2));
    M = STM(:,1+m*(N-1):m*N);

    for i = 1:length(t)
        Phi(:,1+n*(i-1):n*i) = M * ( STM(:,1+m*(i-1):m*i) \ B(:,1+n*(i-1):n*i) );
    end

    % Compute the initial missvector
    b = xf - M * x0;

    % Initial solution (least-squares)
    tic
    x0 = pinv(Phi) * b;
    x0 = reshape(x0, n, []);             % Initial sequence

    % Saturation
    switch (obj.Thruster.q)
        case 'L1'
            q_norm = max(abs(x0), [], 1);

        case 'L2'
            q_norm = sqrt(dot(x0, x0, 1));
            x0 = [x0; ...
                 (-q_norm + obj.Thruster.umin) .* (q_norm < obj.Thruster.umin & q_norm > 0); ...
                 (+q_norm - obj.Thruster.umax) .* (q_norm > obj.Thruster.umax & q_norm > 0)];
            
            if (obj.Thruster.umax == Inf)
                x0 = x0(1:2,:);
            end

            if (obj.Thruster.umin == 0)
                if (size(x0,1) == 3)
                    x0 = [x0(1,:); x0(3:end,:)];
                else
                    x0 = x0(1,:);
                end
            end
    
        case 'Linfty'
            error('Constrained Linfty rendezvous is not yet implemented');
    end

    % Prunning by the Generalized Potter Algoritm 
    
    [sol, obj.Cost] = obj.PVT_pruner(STM, B, x0, obj.Thruster.umax, obj.Thruster.umin, obj.Thruster.p); 
    
    % Polishing 
    if (obj.pol_flag)
        switch (obj.Thruster.q)
            case 'L1'
                ti = max(abs(sol), [], 1) ~= 0;
            case 'L2'
                ti = sqrt(dot(sol, sol, 1)) ~= 0;
            case 'Linfty'
                ti = sum(abs(sol), 1) ~= 0;
        end

        Phi_aux = Phi(:, logical( kron(ti, ones(1,n)) ));
        sol = pinv(Phi_aux) * b;

        dV = zeros(n, N);
        dV(:, logical(ti)) = reshape(sol, n, []);
    else
        dV = reshape(sol, n, []);
    end

    obj.SolveTime = toc;

    % Results
    obj.e(:,1) = b - Phi * reshape(sol, [], 1);          % Final missvector  
    obj.u = dV;                                          % Final rendezvous impulsive sequence
    obj.t = t;                                           % Execution times
    e = obj.e;                                           % Final rendevous error
    u = dV;                                              % Final rendezvous sequence
end