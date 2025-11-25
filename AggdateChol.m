%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 26/11/25
% File: AggdateChol.m 
% Issue: 0 
% Validated:

%% Aggregate-update Cholesky factor %%
% This function updates a Cholesky factor whenever a row is added to the original matrix

function [Lnew] = AggdateChol(L, A, Atilde)
    % Update of the Cholesky factor
    opts.LT = true;

    S12     = A * Atilde.';
    S22     = Atilde * Atilde.';
    Y       = linsolve(L, S12, opts); 
    S22p    = S22 - Y.' * Y; 
    Lnew    = chol(S22p, "lower");
    
    % Final factor
    Lnew = [L zeros(); Y.' Lnew];
end