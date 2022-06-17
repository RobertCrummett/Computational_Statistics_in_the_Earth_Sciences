%% Example 7.19
% Returning to Gossett's wheat data from Example 6.9 and 7.16, use the sign
% rank test to evaluate the null hypothesis that the difference in yeild
% between air- and kiln- dried wheat seeds has a median value of zero.

%% Importing Data
clear; close all; clf
gossett = importdata('gossett.dat');

%% Wilcoxon Signed Rank Test
[p, h, stats] = signrank(gossett(:,1), gossett(:,2));
xihat = stats.signedrank;

% The null hypothesis is accepted. The exact distribution was used for this
% small sample. For bigger samples, an asymptotic normal approximation
% would be used. This can be set by using a keyword in the function too.