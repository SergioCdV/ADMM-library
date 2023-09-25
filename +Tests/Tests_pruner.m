% Sergio Cuevas del Valle
% Date: 25/01/23
% File: PVT_pruner.m 
% Issue: 0 
% Validated: 25/01/23

%% Primer Vector Theory Pruner Tests %%
% This script provides several test cases for the PVT pruner for linear
% systems, as provided by Prussing

clear; 
clc; 

%% Test case I %%
% Time span
t = 0:1e-1:3;               

% Impulse sequence
dV = zeros(2,length(t));
dV(:,1) = [1; 0];
dV(:, t == 1) = [0; 1];
dV(:, t == 2) = [3; 4]; 

B = [1 0; 0 1];

% STM 
Phi = repmat([1 1; 0 1], 1, length(t));
for i = 1:length(t)
    Phi(1,2*i) = t(i);
end

% Pruning 
[dV, cost] = PVT_pruner(Phi, B, dV);

%% Test case II %%
% Time span
t = 0:pi/32:pi/2;               

% Impulse sequence
dV = zeros(1,length(t));
dV(:,1) = -1.5;
dV(:, t == pi/4) = sqrt(2)/2;
dV(:, t == pi/2) = 0.5; 

B = [0; 1];

% STM 
Phi = [cos(t).' sin(t.') -sin(t.') cos(t.')];
Phi = reshape(Phi.', 2, []);

% Pruning 
[dV, cost] = PVT_pruner(Phi, B, dV);
