%% Example 8.12
% Remove the bias from the Weibull fit to earthquake data in Example 7.7

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% Fit Weibull Distribution to Data
% MLE parameters
parmhat = wblfit(quakes);
alpha = parmhat(1);
beta = parmhat(2);

cdf = [sort(quakes) wblcdf(sort(quakes)', parmhat(1), parmhat(2))'];

[h, p, ksstat] = kstest(quakes, cdf);

% The p-value diifers slightly from the one usually used with the bootstrap
% hypothesis test. The sample p-value is 0.9006 and is clearly upward
% biased, although the correction does not change the outcome of the
% hypothesis test as we will see.

n = length(quakes);
b = 999999;
sps = zeros([1 b]);

parfor i = 1:b
    draw = sort(wblrnd(alpha, beta, 1, n));
    cdf = [draw' wblcdf(draw', alpha, beta)];
    [~, ~, ks1] = kstest(draw, cdf);
    sps(i) = ks1;
end

pval = 2*min(sum(sps >= ksstat) + 1, sum(sps < ksstat) + 1)/(b + 1);

% Clearly the sample p-value was much larger than the p-value obtained by
% the K-S bootstrapping. The conclusion is unchanged however.

%% Plotting Monte Carlo Distribution
figure(1);
ksdensity(sps)
hold on
xline(ksstat, 'r--')
hold off
legend('KDE','test statistic','Location','ne')
xlim([0 0.4]); ylim([0 15]);
ylabel('Monte Carlo Distribution')
title 'Monte Carlo Distribution of the Earthquake K-S Statistic'
