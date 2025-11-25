

function [L] = GivensRotation( A, r )
    B = A;                
    [n, m] = size(A);  

    for t = 0:m-r-1
        % Compute Givens parameter 
        k = r + t;
        k = min(k, n);
        b = B(k,k+1);            % l_{r+t, r+t+1}

        if ( b ~= 0 )
            a = B(k,k);           % l_{r+t, r+t}
            denom = hypot(a, b);  % sqrt(a^2 + b^2) estable numéricamente
            c = a / denom;
            s = b / denom;

            for i = k:n
                old1 = B(i,k);
                old2 = B(i,k+1);

                B(i,k)   =  c * old1 + s * old2;
                B(i,k+1) = -s * old1 + c * old2;
            end
        end
    end
    
    % Final output
    L = B;
end