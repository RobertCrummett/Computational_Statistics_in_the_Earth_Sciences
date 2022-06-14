%% Example 6.7
% A meteorologist determines that during the current year there were 80
% major storms in the northern hemisphere. He claims that this represents
% a significant increase over the 100 year average of 70+-2 major storms
% per year.

%% Hypothesis Testing
clear; close all; clf;
% Null hypothesis is that the means are the same. The alternative
% hypothesis is that they are different.

N = 1; % The sample size is one.
mu = 70;
sig = 2;
xbar = 80;

% Z-Score for the observation
zobs = sqrt(N)*(xbar - mu)/sig;

% Critical Value
cv = norminv(0.975);

% Since the abs(zobs) >> cv, this is strong evidence to reject the null
% hypothesis.

% P-value
pval = 2*(1 - normcdf(zobs,0,1));

% Small p-value gives very stong evidence for the alternative hypothesis,
% that this year represents a significant departure from the last 100
% year of data.
