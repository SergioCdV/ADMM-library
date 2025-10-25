%% Modified Chebyshev Picard Iteration Solver %% 
% Sergio Cuevas del Valle
% Date: 19/10/25
% File: MPCI.m 
% Issue: 0 
% Validated: 

%% MCPI solver %%
% This function implements a Modified Chebyshev Picard Iteration TBVP solver

classdef MCPI 
    properties
        N = 1;          % Order of the Chebyshev polynomial approximation
        tol;            % Tolerance to stop the iteration process
        maxIter = 100;  % Maximum number of iterations 
    end

    properties (Access = private)
        Ca;             % Approximation matrix for the dynamics
        Cx;             % Approximation matrix for the state vector
        tau;            % Vector of Chebyshev nodes
    end

    methods 
        % Class constructor
        function [obj] = MCPI(order, tol)
            % Set basic properties 
            obj.N = order + 1; 

            if ( exist('tol', 'var') )
                tol = 2.25E-14;
            end

            obj.tol = tol;

            % Initialize the solver 
            obj = obj.Init();
        end

        % Pre-allocation of matrices
        function [obj] = Init(obj)
            % Constants
            obj.N = obj.N;                              % Number of coefficients in the polynomial expansion
            order = obj.N - 1;                          % Order of the approximation
            Ones = diag(ones(1,obj.N));

            % Weight matrix for the state vector
            W = Ones;                    
            W(1,1) = 1/2;                           

            % Approximation weights for the dynamics
            V = (2 / order) * Ones;              
            V(1,1) = 1/order;                           
            V(obj.N,obj.N) = V(1,1);                           
        
            R = (1/2) ./ ( 1:order );               % Approximation weights for the dynamics
            R = diag([1 R]);                        % Approximation weights for the dynamics

            % Chebyshev polynomials difference
            S = zeros(obj.N);
            S(1,1) = 1; 
            S(1,2) = -1/2;
            S(1,obj.N) = (-1)^(order+1)/(order-1);

            for i = 2:order
                S(1,i+1) = (-1)^(i+1) * 2 / (i^2 - 1);
                S(i,i-1) = +1; 
                S(i,i+1) = -1;
            end
            S(order,end-2)  = +1;
            S(order,end)    = -1;
            S(obj.N, end-1) = +1;

            % Chebyshev coefficients and polynomial matrix
            obj.tau = obj.ClenshawCurtisNodes( order );
            T = obj.ChebyshevPolynomial(order, obj.tau).';

            % Final matrices
            obj.Ca = R * S * (T.' * V);                   % Dynamics approximation matrix
            obj.Cx = T * W;                             % State approximation matrix
        end

        % Compute the physical time corresponding to a given distrbution of CGL nodes
        function [t, J] = PhysicalTime( obj, deltaT, tau )
            % Sanity check
            if ( ~exist('tau', 'var') )
                tau = obj.tau;
            end
            
            % Final output
            t = deltaT * 0.5 * (tau + 1);
            J = 0.5 * ones( 1, length(tau) ); 
        end
    end

    methods (Static)
        % Compute CGL nodes 
        function [tau] = ClenshawCurtisNodes(N)
            i = N:-1:0;
            tau = cos(pi*i/N);
        end
        
        % Compute Chebyshev polynomials
        function [Pn] = ChebyshevPolynomial(order, u)
            % Preallocation of the polynomials 
            Pn = zeros(order+1, length(u)); 
            Pn(1,:) = ones(1, length(u));    % Initialization of the Chebyshev polynomials of the first kind
            Pn(2,:) = u;                     % Initialization of the Chebyshev polynomials of the first kind
    
            % Main computation 
            for i = 2:order
                Pn(i+1,:) = 2 * u .* Pn(i,:) - Pn(i-1,:); % Chebyshev polynomials
            end
        end
    end
end