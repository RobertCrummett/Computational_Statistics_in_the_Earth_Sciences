%% Example 7.5
% Pearson's chi squared test for goodness of fit.

%% Import Data
clear; close all; clf;
bomb = importdata('bomb.dat');

%% Hypothesis Testing
% Observed
o = bomb(:,2);

% Expected
e = poisspdf(bomb(:,1), sum(bomb(:,1).*bomb(:,2))/sum(bomb(:,2)))*sum(bomb(:,2));

% Test Statistic
q = sum((o - e).^2./e);

% P-Value
pval = 1 - chi2cdf(q, length(bomb) - 1 - 1);

% The same result as the multinomial likelihood ratio test.
% This result is again heavily dependant on how the data is binned.
% The problem with this data set is that its bin sizes are too small to
% give adequate sampling, and hence, is not designed for statistical
% testing.