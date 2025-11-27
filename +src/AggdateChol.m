%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 26/11/25
% File: AggdateChol.m 
% Issue: 0 
% Validated:

%% Aggregate-update Cholesky factor %%
% This function updates a Cholesky factor whenever a row is added to the original matrix

function [Lnew, Anew] = AggdateChol(L, A, Atilde)
    % Update of the Cholesky factor
    opts.LT = true;

    S12     = A * Atilde.';
    S22     = Atilde * Atilde.';
    Y       = linsolve(L, S12, opts); 
    S22p    = S22 - Y.' * Y; 

    Anew = [A; Atilde];

    try        
        % Final factor
        Lnew = chol(S22p, "lower");
        O    = zeros( size(L,1), size(Lnew,2) );
        Lnew = [L O; Y.' Lnew];

    catch
        Lnew = chol( Anew * Anew.', "lower");
    end
end