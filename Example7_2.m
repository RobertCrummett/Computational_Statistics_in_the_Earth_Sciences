%% Example 7.2
% Assess whether bomb strikes were random or not.

% If the bomb hits were random, then a Poisson distribution should fit the
% data. The observed occurances are in the second column of the table. The
% expected occurances are given by the Poisson probabilities scales by the
% total number of hits.

%% Importing Data
clear; close all; clf;
bomb = importdata('bomb.dat');

%% Likelihood Ratio Testing for the Multinomial Distribution
% Observed
o = bomb(:,2);
% Expected (for Poisson process)
e = poisspdf(bomb(:,1), sum(bomb(:,1).*bomb(:,2))/sum(bomb(:,2)))*sum(bomb(:,2));

% Test Statistic
loglambda = 2*sum(o.*log(o./e));

% Critical Value
cv = chi2inv(0.975, length(bomb) - 2);

% P-Value
pval = 1 - chi2cdf(loglambda, length(bomb) - 2);

% The null hypothesis that the data is Poisson is accepted strongly.
