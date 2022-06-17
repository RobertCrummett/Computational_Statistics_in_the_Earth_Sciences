%% Example 7.6
% Apply Pearson's chi squared test to the earthquake interval data in
% Example 7.3

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% Goodness of Fit Test
[h, p, stats] = chi2gof(quakes, 'cdf', {@expcdf, expfit(quakes)});

% The null hypothesis is accepted, but the test only has three degrees of
% freedom. These can be increased by writing a MATLAB script like in
% Example 7.3