%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 25/10/25
% File: LorentzProx.m 
% Issue: 0 
% Validated:

%% Lorentz cone proximal operator %%
% This function defines the proximal operator of the Lorentz cone

classdef LorentzProx < src.ProxOperator
    methods
        % Constructor 
        function [obj] = LorentzProx(a, rho)
            if ~exists(rho, 'var')
                rho = 1;
            end

            obj = obj@src.ProxOperator( @(x)LorentzProx.projection(a, x), rho );
        end
    end

    methods (Static)
        % Projection onto the Lorentz cone
        function [y, sigma] = projection(x, t)
            % Points inside the cone
            y = x; 
            sigma = t; 

            % Points outside the cone
            norm_ = vecnorm( x, 2, 1 );
            out_idx = norm_ < -t;
            y(:,out_idx) = zeros(size(y,1),sum(out_idx));
            sigma(out_idx) = 0;

            % Points along the cone
            along_idx = norm_ >= abs( t );
            y(:,along_idx)   = 0.5 * ( 1 + t(along_idx) ./ norm_(along_idx) ) .* x(:,along_idx);
            sigma(along_idx) = 0.5 * ( norm_(along_idx) + t(along_idx) );
        end
    end
end