%% Example 7.18
% Test the null hypothesis (from Example 7.15) that the median of the
% differences between the test and control groups is zero.

%% Importing Data
clear; close all; clf;
remiss = importdata('remiss.dat');
x = remiss(1:21,2);
y = remiss(22:end,2);

%% Wilcoxon Signed Rank Test
[p, h, stats] = signrank(x, y);

% The null hypothesis is rejected at the 0.05 level, as it was for the sign
% test.