%% Example 7.1
% Rutherford, Geiger, and Bateman (1910) present dome of the frist
% statistical data for radioactive decay. Their experiments measured alpha
% decay in polonium over 2068 disjoint 72-second intervals. The columns in
% the Table 7.1 give the number of counts and the number of occurances.

% Based on the physics, the data should be Poisson, and the sample mean
% constitutes the maximum likelihood estimator (MLE) for the parameter in
% the distribution.

%% Data
clear; close all; clf;
rutherford = [57 203 383 525 532 408 273 139 45 27 10 4 1 1];
c = [0 1 2 3 4 5 6 7 8 9 10 11 13 14];

%% Sample Mean
param = sum(rutherford.*c)/sum(rutherford);
n = sum(rutherford);
p = poisspdf(c, param);

% Test Statistic
loglambda = 2*sum(rutherford.*log(rutherford./(n*p)));

% Critical Value
cv = chi2inv(0.975, length(c) - 2);

% P-Value
pval = 1 - chi2cdf(loglambda, length(c) - 2);

% The null hypothesis that the data are Poisson is accepted weakly.
