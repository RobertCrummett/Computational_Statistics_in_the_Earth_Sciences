%% Example 10.1
% For the medical data, suppose that the test data should be 4.1. This is
% very nearly the grand mean of all the data. The mean for each labratory
% can be tested individually against this value with a univariate t test to
% get the p-values, and the hypothesis that all the labratory means
% simultaneously have this value can be evaluated using Hotelling's T2
% test.

%% Importing Data
clear; close all; clf;
tablet = importdata('tablet.dat');

%% Data
[n, m] = size(tablet);
ybar = mean(tablet);
ystd = std(tablet);

%% Univariate t Test
% The test statistic is equal to 4.1
t = sqrt(n)*(ybar - 4.1)./ystd;

pval = 2*(1 - tcdf(abs(t), n - 1));

%% Bonferroni Threshold
A = 0.05;

alpha = A/m;
RejBonIdx = find(pval <= alpha);

%% Benjamini-Hochberg Multiple Hypothesis Test
[RejTh, RejBHIdx] = Benjamini(pval);

%% Multivariate Test Using Hotelling's T2
s = cov(tablet);
t2 = n*(ybar - 4.1*ones(size(ybar)))*inv(s)*(ybar - 4.1*ones(size(ybar)))';

pvalt2 = 2*min(fcdf((n - m)*t2/(m*(n - 1)), m, n - m), 1 - fcdf((n - m)*t2/(m*(n - 1)), m, n - m));

% Since the p-value is below the 0.05 threshold, so the null hypothesis
% that the labratories simultaneously return a value of 4.1 i rejected,
% although not strongly.

%% Confidence Intervals on Sample Means
t2 = (n - 1)*m/(n - m)*finv(1 - A, m, n - m);
ci = zeros(m, 3);

for i = 1:m
    a = zeros(size(ybar));
    a(i) = 1;
    ci(i, :) = [i, a*ybar' - sqrt(t2)*sqrt(a*s*a'/n), a*ybar' + sqrt(t2)*sqrt(a*s*a'/n)];
end

