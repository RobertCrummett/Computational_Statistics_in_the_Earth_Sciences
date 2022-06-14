%% Example 6.6
clear; close all; clf;5
% Popuation mean for p-wave velocity in the mantle is known to be 8km/s
% with a standard deviation of 2km/s.
mu = 8;
sig = 2;

% The sample mean is 7.4km/s. Can it be concluded that the sample mean came
% from a population with mean 8km/s?
xbar = 7.4;
n = 30; % sample size

%% Z-Score Test Statistic
zhat = sqrt(n)*(xbar - mu)/sig;

%% Two Tailed Critical Value and P-Value
cv = norminv(0.975); % two tailed critical value

% a comparison of abs(zhat) and cv leads one to conclude that the null
% hypothesis, that mu = 8km/s, is true.

pval = 2*(1 - normcdf(abs(zhat))); 
% no edividence for alternative hypothesis

%% Power Calculation
muone = 7;
zalow = (sqrt(n)*(xbar - muone)/sig) - cv; % lower hypothesized z-score
zahigh = (sqrt(n)*(xbar - muone)/sig) + cv; % upper hypothesized z-score

% Power of test
power = 1 - normcdf(zahigh) + normcdf(zalow);

% Thus, there is a 19% chance that the null hypothesis is rejected when it
% is false. The power is low because the sample size is small.

% If the power were required to be 0.9, then the sample size would have to
% be approximately 260 for the same experiement with alpha = 0.05