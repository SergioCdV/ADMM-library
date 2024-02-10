%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 09/02/24
% File: Test_LU.m 
% Issue: 0 
% Validated: 

%% Test exercise %% 
% Update an LU through a rank 1 update % 
clear all;

%% Input data
m = 40;
n = 40;
A = rand(m,n);

%% Perturb the matrix 
a = rand(m,1);
p = 2;

%% Compare the LU decompositions
% Original decomposition
[~, ~, P] = lu(A, 'vector');
A = A(P,:);
[L, U, ~] = lu(A);

A(:,p) = a;

% Update 
[L2, U2] = Solvers.FT_update(L, U, a, p);

%% Check the nullspace 
% Generate a rank-defficient matrix
m = 2;
n = 3;
A = rand(m, 2*n+1); 

% Direct nullspace through SVD
N = null(A);

% Nullspace from Lu
[~, ~, P] = lu(A, 'vector');
A = A(P,:);
[L, U, ~] = lu(A);

m = size(U,1);
d = size(A,2) - m;
N2 = [-U(:,1:m) \ U(:,m+1:end); eye(d)];

x = ones(size(N,2),1); 

norm( A * N2 * x)

%% Solve linear system 
% Problem data
A = rand(m,m);
b = rand(m,1);

% LU factorization
[~, ~, P] = lu(A, 'vector');
A = A(P,:);
[L, U, ~] = lu(A);

% Solve the system
x = LUsolve(L,U,b);

% Results
norm(b - A * x);

%% Auxiliary functions 
% Solve a linear system using LU
function [x] = LUsolve(L,U,b)
    y = FTRAN(L,b);
    x = BTRAN(U,y);
end

% Forward substitution
function [x] = FTRAN(L, b)

    x = zeros(size(L,1),1);
    for i = 1:size(L,1)
        x(i,1) = ( b(i,1) - dot(L(i,1:i-1), x(1:i-1,1)) ) / L(i,i);
    end

end

function [x] = BTRAN(U, b)
    
    n = size(U,1);
    x = zeros(size(U,1),1);
    for i = n:-1:1
        x(i,1) = ( b(i,1) - dot(U(i,i+1:n), x(i+1:n,1)) ) / U(i,i);
    end

end