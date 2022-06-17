%% Example 6.28
% Apply the Rayleigh test to the wind direction data.

%% Rayleigh Test
clear; close all; clf;
wind = importdata('wind.dat');

% Rose Diagram
figure(1);
polarhistogram(wind, 'NumBins', 25);
title 'Rose diagram of directional wind data from 2004 - 2005'

s = sum(sind(wind));
c = sum(cosd(wind));
rbar = sqrt(s^2 + c^2)/length(wind);

% Score Statistic
rayl = 2*length(wind)*rbar^2;

% P-Value
pval = 1 - chi2cdf(rayl, 2);

% The null hypothesis that the data are uniformly distributed is rejected,
% which is no big suprise veiwing the figure.