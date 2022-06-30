%% Example 10.8
% Compute the canonical correlation among the variables.

% Unlike PCA, canonical correlation focus is on the correlation between two
% sets of random variables rather than the correlation within a single data
% set.

%% Importing Data
clear; clf; close all;
sake = importdata('sake.dat');

%% Data
y = sake(:, 1:2);
x = sake(:, 3:10);

%% Canonical Correlation
[a, b, r, u, v, stats] = canoncorr(x, y);

% The canonical correlations (r) are moderate, but they are not significant, as
% can be shown by examining stats:

disp(stats.pF) % p-value for F approximation to Wilk's lambda
disp(stats.pChisq) % p-value for chi squared approximation

% In both approximations, the null hypothesis that the data are
% uncorrelated is accepted. There is no discernable correlation in this
% data.