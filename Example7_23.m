%% Example 7.23
% Compute the Spearman rank correlation test on the alcohol tabacco data
% from Example 6.17

%% Importing Data
clear; close all; clf;
data = importdata("alctobacc.dat");

%% Spearman Rank Correlation Test
[rhat, p] = corr(data(:,1),data(:,2),'type','Spearman');

% The null hypothesis that the data are uncorrelated is accepted at the
% 0.05 level. Compared to Example 6.17, the p-value has fallen and the
% correlation has risen. This indicates that the presence of the outlier is
% still felt, but the test is less sensitive to it.

% Once the outlier (the last data point) is removed, the alternative
% hypothesis is strongly accepted - similar to Example 6.17.

%% Spearman Rank Correlation Test w/o Outlier
% Data has been 'cleaned'
[rhat_c, p_c] = corr(data(1:end-1,1),data(1:end-1,2),'type','Spearman');