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

%% Pertub the matrix 
a = rand(m,1);
p = 2;

%% Compare the SVD decompositions
% Original decomposition
[~, ~, P] = lu(A, 'vector');
A = A(P,:);
[L, U, ~] = lu(A);

A(:,p) = a;

% Update 
[L2, U2] = Solvers.FT_update(L, U, a, p);
