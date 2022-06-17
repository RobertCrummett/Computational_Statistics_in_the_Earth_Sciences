%% Example 7.4
% Pearson's chi squared test for goodness of fit.

%% Data
clear; close all; clf;
rutherford = [57 203 383 525 532 408 273 139 45 27 10 4 1 1];
c = [0 1 2 3 4 5 6 7 8 9 10 11 13 14];

%% Hypothesis Testing
param = sum(rutherford.*c)/sum(rutherford);
n = sum(rutherford);
p = poisspdf(c, param);

% Test Statistic
q = sum((rutherford - n*p).^2./(n*p));

% P-Value
pval = 1 - chi2cdf(q, length(c) - 1 - 1);

% The null hypothesis that the data are Poisson is rejected.

%% Examine the effect of class boundaries...
c = [c(1:10) 10];
rutherford = [rutherford(1:10) 16];

param = sum(rutherford.*c)/sum(rutherford);

n = sum(rutherford);
p = poisspdf(c, param);

% Test Statistic
q2 = sum((rutherford - n*p).^2./(n*p));

% P-Value
pval2 = 1 - chi2cdf(q2, length(c) - 1 - 1);

% The null hypothesis that the data are Poisson is accepted.