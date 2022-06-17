%% Example 7.17
% Use the Wilcoxon signed rank test to evaluate the null hypothesis that
% the number of days per year that an oceanographer spends at sea is 30.

% A downside of the sign test is that only the sign of the difference
% between a sample and the postulated median enters it. The Wilcoxon signed
% rank test also incorporates the avsolute values of the differences.

%% Data
clear; close all; clf;
x = [54 55 60 42 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];

%% Wilcoxon Signed Rank Test
mustar = 30;
[p, h, stats] = signrank(x, mustar);

% The null hypothesis is weakly rejected at the 0.05 level, as was the case
% for the t test in Example 6.8, but not for the sign test in Example 7.14.