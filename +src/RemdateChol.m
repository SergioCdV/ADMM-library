

% Update Cholesky factor whenever a row is deleted
function [Lnew] = RemdateChol(L, k)
    % Make a copy of the previous matrix
    H = L;

    % Sort rows to ease the process
    r = sort(k);

    for i = 1:length(r)
        % Delete the r-th row of the original matrix 
        idx = r(i) - (i-1);
        H(idx,:) = [];
    
        % Do Givens rotations along the required non-diagonal terms 
        H = GivensRotation( H, idx );

        % Finally, delete the last row and column 
        H = H(:,1:end-1);
    end
    
    % Final output
    Lnew = H;
end