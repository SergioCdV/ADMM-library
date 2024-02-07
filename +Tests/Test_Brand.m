%% ADMM Library %% 
% Sergio Cuevas del Valle
% Date: 06/02/24
% File: Test_Brand.m 
% Issue: 0 
% Validated: 

%% Test exercise %% 
% Update an SVD through Brand's method % 
clear all;

%% Input data
m = 2;
n = 2;
A = rand(m,n);

%% Compute the original SVD 
[U,S,V] = svd(A);

%% Pertub the matrix 
a = rand(m,1);
b = [zeros(n-1,1); 1]; 

A = A + a*b.';

%% Compare the SVD decompositions
% Direct SVD
A = [A a];
[S1,V1,D1] = svds(A); 

% Updated SVD 
[S2,V2,D2] = Solvers.Brand_SVD(U, S, [V; zeros(1,size(V,2))], a, [zeros(1,size(V,2)) 1].');

%% Results
Ar = S2*V2*D2.';
error = norm(Ar-A,"fro");
norm(D2 * D2.')
norm(S2 * S2.')





