%% Example 8.10
% From the cloudseeding data of Example 6.12, compute a bootstrap test for
% the null hypothesis that there is no effect from cloud seeding against
% the alternative that there is.

%% Importing Data
clear; close all; clf;
clouds = importdata('cloudseed.dat');
rng default

%% Two Sample Test for the Equality of Means
b = [99 999 9999 99999];

% Two Sample t Statistic
that = (mean(clouds(:,1)) - mean(clouds(:,2)))./sqrt(var(clouds(:,1)) + var(clouds(:,2)));

% Grand Mean of the Data
gmean = (mean(clouds(:,1)) + mean(clouds(:,2)))/2;

%% Bootstrapping
for i = 1:length(b)
    bmat1 = bootstrp(b(i), @(x) [mean(x) var(x,1)], clouds(:,1) - mean(clouds(:,1)) + gmean);
    bmat2 = bootstrp(b(i), @(x) [mean(x) var(x,1)], clouds(:,2) - mean(clouds(:,2)) + gmean);
    boot = (bmat1(:,1) - bmat2(:,1))./sqrt(bmat1(:,2) + bmat2(:,2));
    pval(i) = 2*min(sum(boot > that),sum(boot <= that))./b(i);
end

% The null hypothesis is rejected, and there is a difference in rainfall
% between seeded and unseeded clouds. This is in agreement with the rank
% sum test but in contrast with the two sample t test.

%% Plotting Bootstrap Data
figure(1);
ksdensity(boot)
xline(that,'r--')
title 'KDE of the bootstrap distirbution of cloud data resampled 99999 times'
legend('KDE','test statistic','Location','ne')

% The distribution is asymmetric due to the logarithmic nature of the data.
% The sample test statistic is -0.3919. A nearly identical result is
% obtained after taking logs of the data. Thus, this illustrates how the
% bootstrapping method is robust against departures from Gaussianity.
