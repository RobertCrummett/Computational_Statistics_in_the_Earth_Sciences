%% Example 6.18
% There are ten representative measuements per labratory. The goal is to
% determine whether there are systematic differences between the labratory
% results.

%% Import Data
clear; close all; clf;
tablet = importdata('tablet.dat');

%% ANOVA, analysis of variance
[p, anovatab] = anova1(tablet);

N = size(tablet, 1);
M = size(tablet, 2);

% test statistic
F = anovatab{2,5};

% P-Value
pval = 2*min(1 - fcdf(F, M - 1, M*(N - 1)), fcdf(F, M - 1, M*(N - 1)));

% The null hypothesis is rejected. The means are not all the same.