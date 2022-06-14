%% Example 6.17
% The correlation coefficent allows two populations drawn from a bivariate
% normal distribution to be compared. The null hypothesis under the test
% says that there is no correlation between the data.

% Assess whether there is significant correlation between these two data
% sets and investigate whether the result is unduly influenced by a single
% datum.

%% Import Data
clear; close all; clf;
alctobacc = importdata('alctobacc.dat');

%% Remove Outlier
% To remove the outlier, uncomment line 17.

%alctobacc = alctobacc(1:end-1,:);

%% Scatter Plot of Data
scatter(alctobacc(:,2),alctobacc(:,1),'r+')
xlabel 'Tobacco Spending'; ylabel 'Alcohol Spending'; 
title 'British Household Spending in 11 Regions of the United Kingdom';
ylim([3.5 6.5])

%% Correlation Coefficent
X = alctobacc(:,1);
Y = alctobacc(:,2);
N = length(alctobacc(:,1));
Xbar = mean(X);
Ybar = mean(Y);

% Sample Correlation Coefficent
rhat = (sum((X - Xbar).*(Y - Ybar)))/sqrt(sum((X - Xbar).^2)*sum((Y - Ybar).^2));

Corr = 1/2 + gamma((N-1)/2)*rhat*Hypergeometric2f1(1/2,2-(N/2),3/2,rhat^2)/(sqrt(pi)*gamma(N/2 - 1));

% P-Value of Data
pval = 2*min(Corr, 1- Corr);

% The high p-value is no evidence for the alternative hypothesis.
% When the data is trimmed, notice the changes!

% To do this correctly, after deleting the outlier, the null distribution
% should be taken as the truncated form
