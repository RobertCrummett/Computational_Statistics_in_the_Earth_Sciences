%% Example 7.14
% Preform the sign test for the median (rather than just the mean) number
% of days an oceanographer spends at sea per year.

%% Data
clear; close all; clf;
x = [54 55 60 42 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];

%% Sign Test
mustar = 30;
[p, h, stats] = signtest(x, mustar);
xihat = stats.sign;

% The null hypothesis that the median number of days an oceanographer
% spends at sea per year is 30 is accepted.