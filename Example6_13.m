%% Example 6.13 & 14 - Chi Squared Testing
% The chi quared test applies when the distribution for the test statistic
% follows the chi squared distirbution under the null hypothesis. In its
% simplest form, the chi squared test is used to determine whether the
% variance computed from a sample has a specified value when the underlying
% distirbution of the sample is Gaussian.

%% Importing Data
clear; close all; clf;
data = importdata('battery.dat');

variance = 5;
sig2 = var(data, 1); % biased sample variance
N = length(data); % sample size

%% Hypothesis Testing
% Test Statistic
chi2hat = N*sig2/variance;

% Critical Value
cv = chi2inv(0.975, N - 1); % Two-tailed test, alpha = 0.05

% The test statistic is less than the critical value, so the null
% hypothesis is accepted.

pval = 2*(1 - chi2cdf(chi2hat, N - 1));

% There is weak evidence for the alternative hypothesis.

% Recomputing for the upper tail test, it is found that the null hypothesis
% fails with strong evidence for the alternative hypothesis!

%% Plotting Power Curve
% The observed variance ratio
varrat_obs = sig2/variance;

% Two-tailed test boundaries
xlo = chi2inv(0.025, N - 1);
xhi = chi2inv(0.975, N - 1);

% Upper tailed test boundary
xhi1 = chi2inv(0.95, N - 1);

% Power Calculation for Both Tests
varrat = 1:0.01:2;
delta = N*(varrat - 1); % Non-centrality parameter
power = 1 - ncx2cdf(xhi, N - 1, delta) + ncx2cdf(xlo, N - 1, delta);
power1 = 1 - ncx2cdf(xhi1, N - 1, delta);

% Plotting
figure(1);
plot(varrat, power)
hold on
plot(varrat, power1)
hold off

xlabel 'Variance Ratio';
ylabel 'Power';
title 'Power Curves for Two-Sided and Upper-Tail Test on Battery Data';
legend('Two-Tailed Test','Upper-Tail Test',Location='nw')