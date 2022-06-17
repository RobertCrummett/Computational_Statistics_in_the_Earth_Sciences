%% Example 7.21
% From the cloud seeding data, compute the p-value for the rank sum test to
% assess the null hypothesis that the medians of the rainfall from
% unseeded and seeded clouds are the same.

%% Importing Data
clear; close all; clf;
clouds = importdata('cloudseed.dat');

%% Rank Sum Test
[p, h] = ranksum(clouds(:,1), clouds(:,2));

% The null hypothesis is rejected, unlike the two-sample t test in Example
% 6.12.

%% Rank Sum Test on Logs of the Data
[plog, hlog] = ranksum(log(clouds(:,1)), log(clouds(:,2)));

% The p-value is unchanged becaause the ranks are not affected by monotone
% tranformations like the logarithm.

% This is an example of the robustness of the rank sum test to departures
% from normality.