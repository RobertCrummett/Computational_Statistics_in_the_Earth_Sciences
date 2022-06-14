%% Example 6.8
% The t test applies when the distribution for the test statistic follows
% Student's t distribution under the null hypothesis. It requires that the 

%% Hypothesis Testing
% Null hypothesis is that mu = 30 days. The alternative hypothesis is that
% the mean is not 30 days at sea.
clear; close all; clf;

mu = 30;
xhat = [54 55 60 42 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];
N = length(xhat);
xbar = mean(xhat); % sample mean
sig = std(xhat,1); % sample baised standard deviation

% Test statistic for baised standard deviation
that = sqrt(N - 1)*(xbar - mu)/sig;

% Two-Tailed Critical Value
cv = tinv(0.975, N - 1);

% Since the test statistic is greater than the critical value, the null
% hypthesis is rejected.

% P -Value
pval = 2 * (1 - tcdf(that, N - 1));

% Small p-value is strong evidence for the alternative hypothesis.

%% Plotting Power Curves of the One-Sample t Test
mualt = 10:0.1:50;
effect = (mualt - mu)/std(xhat);
% The population standard deviation in the effect has had to be replaced by
% the sample standard deviation.

% Non-centrality parameter
delta = sqrt(N)*effect;
power = 1 - nctcdf(that, N - 1, delta) + nctcdf(-that, N - 1, delta);

figure(1);
plot(effect, power, 'k', 'LineWidth', 2)

xlabel 'Effect Size'; 
ylabel 'Power'; 
title 'Power Curve for Single Sample t Test';
ylim([0 1])

% Effect size is a measure of the difference between the population and
% postulated means normalized by the standard error.

% The minimum power is 0.05 at xbar = 30.

% A test whose power is no smaller than the significance level us unbiased,
% and in fact, the t test is a uniformly most powerful unbiased test.