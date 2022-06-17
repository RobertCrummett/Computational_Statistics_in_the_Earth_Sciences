%% Example 7.12
% Compute A-D test on the random normal data from Example 7.8, after it has
% been contaminated.

%% Data
clear; close all; clf;
rng default;
data = normrnd(0, 1, 1000, 1);

% Contamination
data([1000, 999, 998, 997]) = [6.5 6.0 5.5 5];

%% A-D Test
[h, p, adstat, cv] = adtest(data);

% The test strongly supports the alternative hypothesis that the data are
% not normal. This illustrates the sensitivity of the Anderson-Darling in 
% the tails as compared to the Kolmogorov-Smirnov test.