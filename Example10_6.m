%% Example 10.6
% Analyze Siotani et al.'s (1963) data set containing two response
% variables and eight predictor variables. The predictive power of the data
% will be evaluated.

%% Importing Data
clear; clf; close all;
sake = importdata('sake.dat');

%% Data
y = sake(:, 1:2);
x = [sake(:, 3:10) ones(length(y), 1)];
[n, m] = size(x);

%% Multivariate Regression
beta = x\y;

e = (y - x*beta)'*y;
h = beta'*x'*y - n*mean(y)'*mean(y);

[pval, lambda] = WilksLambda(h, e, 2, m - 1, n - m - 2);

% The p-value is large. The null hypothesis that all of the regression
% coefficents are zero is accepted. There is no point in doing any further
% analysis.