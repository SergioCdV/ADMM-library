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