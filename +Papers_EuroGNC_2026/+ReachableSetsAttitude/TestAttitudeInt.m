%% Optimal Linear Slew via ADMM %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: TestAttitudeInt.m 
% Issue: 0 
% Validated: 

%% Test: attitude integration %% 
% Test the MPCI solver on an attitude integration test %

%% Initial data
% Initial conditions 
s0 = [0 0 0 1 1 0 1].';                     % Initial quaternion and initial angular velocity 
sf = [0 sqrt(2)/2 0 sqrt(2)/2 1 0 1].';     % Final conditions

I = eye(3);                 % Inertia matrix of the body 
Iinv = inv(I);              % Inverse of the inertia matrix

t0 = 0;                     % Initial integration time
tf = 2*pi;                  % Final integration time

%% Integration 
% Time span 
delta_t = tf - t0; 

% MCPI integration 
Solver = MCPI( 40, 1E-10 );
tspan = Solver.PhysicalTime( delta_t );

AttitudeSystem = Systems.System( @(t,s)dynamics( I, Iinv, t, s, zeros(3,1) ), Systems.ProblemTypes.TBVP, [s0 sf]);

% Solve the problem
x0 = repmat(s0, 1, Solver.N);
[t, Smcpi] = Solver.Solve( AttitudeSystem, x0, delta_t );

% ODE45 integration 
options = odeset('RelTol', 2.25E-14, 'AbsTol', 1E-22);

[t, Sode] = ode45(@(t,s)dynamics(I, Iinv, t, s, zeros(3,1)), tspan, s0, options);

%% Comparison of results 
error = Smcpi(:, 5:7) - Sode(:,5:7); 

figure
plot( t, vecnorm(error, 2, 2) )
xlabel('t')
ylabel('error')

figure 
plot(t, Smcpi)
grid on;

%% Auxiliary functions 
% Attitude dynamics problem
function [ds] = dynamics(I, Iinv, t, s, u)
    % State variables
    q = s(1:4,:);           % Attitude quaternion
    omega = s(5:7,:);       % Angular velocity, body frame, rad/s

    % Attitude dynamics 
    H = I * omega;
    ds(5:7,:) = Iinv * ( cross(H, omega) + u );

    % Attitude kinematics 
    Omega = [omega; zeros(1,size(omega,2))];            % Angular velocity as a pure quaternion 
    ds(1:4,:) = 0.5 * quaternion_product( Omega, q );
end

% Quaternion product in vectorized form 
function [output] = quaternion_product( left_q, right_q )
    % Pre-allocation 
    m = size(right_q,2);
    n = size(left_q,2); 
    output = zeros(4, n);

    for i = 1:n
        % Left isoclinic
        Sigma = [ [left_q(4,i) * eye(3) + hat_map(left_q(1:3,i)), left_q(1:3,i)]; -left_q(:,i).'];
        output(:,i) = Sigma * right_q(:,i);
    end
end

% Hat map 
function [S] = hat_map(omega)
    S = zeros(3,3); 

    S(1,2) = -omega(3); 
    S(1,3) = +omega(2); 
    S(2,3) = -omega(1);

    S(2,1) = -S(1,2);
    S(3,1) = -S(1,3);
    S(3,2) = -S(2,3);
end