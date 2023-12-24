%% Optimal rendezvous by ADMM %% 
% Sergio Cuevas del Valle
% Date: 24/12/23
% File: ConeProj.m 
% Issue: 0 
% Validated: 

%% Cone Projection %% 
% Projection onto the problem's cone

% Inputs:  - vector indices, determining the size of each primer vector 
%          - q, the thruster maximum control authority norm
%          - b, the linear cost function (miss-vector)
%          - z, the optimization variables
% Outputs: - vector z, the update impulsive sequence

function [cone] = ConeCheck(indices, q, b, z)
    % Primer vector update
    lambda = z(1);     % Lagrange multiplier
    p = z(2:end);      % Primer vector

    % Lagrange multiplier update (projection on a half space)
    cone = zeros(length(z),1);

    cone(1) = lambda > 0;
    start_ind = 1;
    for i = 1:length(indices)
        sel = start_ind:indices(i);

        % Projection of the primer vector on the unit lq-ball
        switch (q)
            case 'L1'
                cone(1+i) = sum(abs(p(sel))) > 1;
            case 'L2'
                cone(1+i) = norm(p(sel)) > 1;
            case 'Linfty'
                cone(1+i) = max(abs(p(sel))) > 1;
        end
                   
        start_ind = indices(i) + 1;
    end
end