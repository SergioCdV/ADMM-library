%% Optimal control by ADMM %% 
% Sergio Cuevas del Valle
% Date: 27/10/25
% File: InputReconstruction.m 
% Issue: 0 
% Validated: 

%% Neustadt Input Reconstruction Algorithm %% 
% Recover the action sequence from the primer vector %

% Inputs:  - vector t, the 1 x n array of impulsive epochs considered
%          - array b , the m x 1 missvector of initial conditions
%          - array Y, the m x N * n array of STM convoluted with the input matrix of the system
%          - array p, the n x N array of computed primer vector
%          - scalar epsilon, numerical tolerance to detect the impulses epochs

% Outputs: - vector t, of dimensions 1 x N, at which the control is to be applied (maneuver execution times)
%          - array u, of dimensions n x N, the control law to be applied (maneuver magnitudes)

function [t, u] = ImpulseReconstruction(t, b, Y, p_norm, epsilon)
    % Constants
    N = size(t,2);        % Number of control epochs
    n = size(Y,1) / N;    % Control dimension

    % Impulsive epochs
    imp_opp = abs(p_norm - 1) <= epsilon;

    if ( sum(imp_opp) > 1E3 )
 % TODO: analysis based on the derivative of the primer vector
%         dp = p_norm - 1;
%         d_dp  = gradient(dp); 
%         dd_dp = gradient(d_dp);
%         dp = abs(dp);
%         [~, pos] = sort(dp);
%         index = zeros(1,length(p_norm));
% 
%         k = 0;
%         for i = 1:length(pos)
%             if (t(pos(i)) == t(1) && sign(d_dp(1)) < 0)
%                 index(pos(i)) = 1;
%                 k = k+1;
%             elseif (t(pos(i)) == t(end) && sign(d_dp(end)) > 0)
%                 index(pos(i)) = 1;
%                 k = k+1;
%             elseif (sign(dd_dp(pos(i))) < 0)
%                 index(pos(i)) = 1;
%                 k = k+1;
%             end
% 
%             if (k == m)
%                 break;
%             end
%         end
%         
%         index = logical(index);
%         t_pruned = t_pruned(index);
%         index = kron(index, ones(1,n));

    elseif ( ~isempty(imp_opp) )
        % Check actuation epochs
        Ones = ones(1,n);                   % Vectors of 1
        t     = t( logical(imp_opp) );      % Pruned impulsive epochs
        index = kron( imp_opp, Ones );      % Indices of the epochs
        index = logical( index );           % Indices of the epochs

        % Compute the maneuver sequence via solving the corresponding linear system
        dv = Y(index,:).' \ b; 
        u = reshape( dv, n, [] );   

    else
        % Null actuation
        u = zeros(n, N);
    end
end