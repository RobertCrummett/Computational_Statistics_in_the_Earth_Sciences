%% Example 7.16
% Use the sign test on Gossett's data to evaluate the null hypothesis that
% the difference in yeild between air- and kiln-dried wheat seeds has a
% median value of zero.

%% Importing Data
clear; close all; clf
gossett = importdata('gossett.dat');

%% Sign Test
[p, h, stats] = signtest(gossett(:,1), gossett(:,2));
xihat = stats.sign;

% The null hypothesis is accepted at the 0.05 level, which also pertains to
% the t test in Example 6.9