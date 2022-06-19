%% Example 7.25
% Non-parametric ANOVA, the Kruskal-Wallis test. Use it on the labratory
% results from Example 6.17.

% From Kvam & Vidakovic (2007, page 108).

%% Importing Data
clear; close all; clf;
tablet = importdata('tablet.dat');

%% Kurskal-Wallis Test
[p, anovatab] = kruskalwallis(tablet);

% The null hypothesis that the medians are the same is rejected at the 0.05
% level. This is the same result as obtained with ANOVA in Example 6.17.