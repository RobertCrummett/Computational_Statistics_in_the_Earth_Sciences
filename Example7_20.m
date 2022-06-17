%% Example 7.20
% Compute the rank sum test on the Michelson speed of light data from
% Example 6.10.

%% Importing Data
clear; close all; clf;
michelson1 = importdata('michelson1.dat');
michelson2 = importdata('michelson2.dat');

%% Rank Sum Test
[p, h, stats] = ranksum(michelson1, michelson2);

% The null hypothesis that the medians of the two populations are the same
% is rejected. Either the speed of light changed or there was systematic
% error in the meaasurements.

% The normal approximation was used to compute the p-value above. The exact
% p-value can be obtained, although it takes longer to calculate.

%% Exact Rank Sum Test
[p_exact, h_exact, stats_exact] = ranksum(michelson1, michelson2, 'method', 'exact');