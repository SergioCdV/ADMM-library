

function [L, U] = FT_update(L0, U0, b, p)
    % Get the permutation matrix
    P1 = eye(size(L0,1));
    L1 = P1;
    P1(:,p:end) = P1(:, [p+1:size(U0,2) p]);
    P1 = P1.';

    R = U0;
    R(1,p) = b(1,1) / L0(1,1);
    for i = 2:size(R,1)
        R(i,p) = ( b(i,1) - dot(L0(i,1:i-1), R(1:i-1,p)) ) / L0(i,i);
    end

    R = (P1 * R) * P1.';
    U1 = R;

    % Compute the small update matrices
    for i = 2:size(L1,2)-1
        L1(end,i) = ( R(end,i) - dot(U1(1:i-1,i), L1(end,1:i-1)) ) / U1(i,i);
    end
    
    U1(end,:) = zeros(1,size(U1,2));
    U1(end,end) = R(end,end) - dot(L1(end,1:end-1), U1(1:end-1,end));

    % Update
    L = L0 * P1.' *L1;
    U = U1 * P1;
end