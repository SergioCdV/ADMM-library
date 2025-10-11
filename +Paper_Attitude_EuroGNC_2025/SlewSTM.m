%% Optimal Linear Slew via ADMM %% 
% Sergio Cuevas del Valle
% Date: 08/10/25
% File: SlewSTM.m 
% Issue: 0 
% Validated: 

%% Slew STM %% 
% Compute the Euler-Poinsot STM in the axisymmetric problem %

function [Phi, PhiInv] = SlewSTM(b, omega0, r0, delta_t)
    % Pre-allocation 
    QuatPhi = zeros(4,4);

    % STM for the attitude kinematics
    Omega = 0.5 * [omega0 + b * r0, omega0 - b * r0];                                  % Fundamental frequencies of motion
    delta_phi = Omega * delta_t;

    cos_omega = cos(delta_phi);                                                        % Pre-allocation of trigonometric functions
    sin_omega = sin(delta_phi);

    QuatPhi(1:2,1:2) = [cos_omega(1) -sin_omega(1); sin_omega(1) cos_omega(1)];        % Motion due to the first frequency
    QuatPhi(3:4,3:4) = [cos_omega(2) -sin_omega(2); sin_omega(2) cos_omega(2)];        % Motion due to the second frequency
        
    % STM for the angular velocity
    delta_phi = b * r0 * delta_t;

    cos_bt = cos(delta_phi);
    sin_bt = sin(delta_phi);
    OmegaPhi = [cos_bt -sin_bt 0; sin_bt cos_bt 0; 0 0 1];

    % Final STM (block-wise constructed)
    Phi = [QuatPhi zeros(4,3); zeros(3,4) OmegaPhi];            % STM of the complete state vector
    PhiInv = [QuatPhi.' zeros(4,3); zeros(3,4) OmegaPhi.'];     % Inverse of the STM
end