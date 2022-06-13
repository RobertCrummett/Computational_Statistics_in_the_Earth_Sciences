%% Example 5.29
% Equations can be solved numerically, or just numerically approximated.
% MLE for binomial distribution phat = 0.25, for j = 5 and N = 20.

clear; close all; clf;

lowerBound = binoinv(0.025, 20, 0.25);
upperBound = binoinv(0.975, 20, 0.25);

% Not symmetric about 5 because the binomial distribution is skewed when p
% does not equal 0.5.

% These answers are approximate because binoinv() returns the least integer
% such that the binomial cdf is evaluated at x equals or exceeds the
% probabilities 0.025 or 0.975. They are conservative confidence intervals.