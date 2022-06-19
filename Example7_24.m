%% Example 7.24
% Compute Kendall's Tau for the alcohol and tabacco consumption data from
% Example 6.17.

%% Importing Data
clear; close all; clf;
data = importdata('alctobacc.dat');

%% Kendall's Tau
[rho, p] = corr(data(:,1),data(:,2),'type','Kendall');

% This is comparable to the result using the Spearman statistic in Example
% 7.23. The uncleaned data still results in zero correlation. After
% removing the data point, the null hypothesis is rejected once again.

%% Kendall's Tau w/o Outlier
[rho_c, p_c] = corr(data(1:end-1,1),data(1:end-1,2),'type','Kendall');

% Strong rejection of the null hypothesis.