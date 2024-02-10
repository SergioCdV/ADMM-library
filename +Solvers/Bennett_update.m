

function [L, U, stable_flag] = Bennett_update(L, U, u, v)
    % Constants 
    m = size(U,1);
    M = size(U,2);
    stable_flag = true;

    % Row-wise algorithm for efficiency
    for i = 1:m
        for j = 1:i-1
            u(i) = u(i) - u(j)*L(i,j);
            L(i,j) = L(i,j) + v(j)*u(i);
        end

        U(i,i) = U(i,i) + u(i)*v(i);
        v(i) = v(i) / U(i,i);

        for j = i+1:M
            U(i,j) = U(i,j) + v(j)*u(i);
            v(j) = v(j) - v(i)*U(i,j);
        end

        % Check for stability
        if (abs(U(i,i)) < max(abs(U(i,2))))
            stable_flag = false;
            break; 
        end
    end
end