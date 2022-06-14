%% Example 6.9
% The t test can be used to compare two paired samples. This is frequently
% used in designed experiments in which a population is sorted into groups
% according to some characteristic.

% In his seminal 1908 paper, W. S. Gossett (aka Student) presented the
% paired results of yeild from test plots of wheat, where half the plots
% were planted with air-dried seed and the remaining plots were okanted
% with kiln-dried seeds. The paired t test can be used to compare the two
% populations for a common mean.

%% Importing Data
clear; close all; clf;
gossett = importdata("gossett.dat");

%% Paired t Test
[h, p, ci, stats] = ttest(gossett(:,1), gossett(:,2));

% The null hypothesis is accepted (h = 0) at alpha = 0.05.
% The probability value is higher than 0.1, which is no evidence to accept
% the alternative hypothesis. The means are the same.