%% Example 6.26
% The file contains number of tornados in a single stat from 1959 to 1990.
% Consider the null hypothesis that the data are Poisson against the
% alternative that they are too dispersed to be Poisson.

% This is called the Poisson Dispersion Test.

%% Import Data
clear; close all; clf;
x = importdata('tornados.dat');

figure(1);
plot(x(:,1),x(:,2),'*k')
xlabel 'Year'; ylabel 'Number of Tornados';
title 'Number of Tornadoes as a Function of Time'

%% Likelihood Testing
% The null hypothesis is that the data are distributed according to a
% Poisson process. The alternative is that they all have different
% dispersion parameters.

torn = x(:,2);
N = length(torn);

% Test Statistic for Poisson Dispersion Test
test = 2*sum(torn.*log(torn/mean(torn)));

% Critical Value
cv = chi2inv(0.975, N - 1);

% P-Value
pval = 1 - chi2cdf(test, N - 1);

% The null hypothesis is accepted, and the data appear to be Poisson.