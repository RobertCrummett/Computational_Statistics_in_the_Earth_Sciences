%% Example 6.12
% In 1950s and 1960s cloud seeding data was collected over several
% countries. This experiment randomly selected 26 clouds for seeding and 26
% clouds for not seeding. Then they measured the acre-feet of rain from all of
% the clouds.

%% Imprting Data
clear; close all; clf;
clouds = importdata('cloudseed.dat');

%% Cloud Seeding q-q Plots

% For the original cloud seeding data
figure(1);
subplot(1,2,1);
qqplot(clouds(:,1))
axis square
title 'Unseeded Cloud Data'
subplot(1,2,2);
qqplot(clouds(:,2))
axis square
title 'Seeded Cloud Data'

% For the logs of cloud seeding data
figure(2);
subplot(1,2,1);
qqplot(log(clouds(:,1)))
axis square
title 'Log Unseeded Cloud Data'
subplot(1,2,2);
qqplot(log(clouds(:,2)))
axis square
title 'Log Seeded Cloud Data'

%% Two-Sample t Test
% Since the q-q plots show departures from Gaussianity, there is reason to
% question whether the t test will preform well or not.

[h, p] = ttest2(clouds(:,1),clouds(:,2));

% The null hypothesis that the clouds produce different amounts of rain
% is weakly accepted.

%% Hypothesis Testing Continued
xbar1 = mean(clouds(:,1));
xbar2 = mean(clouds(:,2));
n1 = length(clouds(:,1));
n2 = length(clouds(:,2));
sig1 = std(clouds(:,1),1);
sig2 = std(clouds(:,2),1);

% Test Statistic
that = (xbar1 - xbar2)/sqrt(sig1^2/(n1 -1) + sig2^2/(n2 -1));
nu = ((sig1^2/(n1 - 1) + sig2^2/(n2 - 1))^2)/((sig1^2/(n1 - 1))^2/n1 + (sig2^2/(n2 - 1))^2/n2) - 2;

% Critical Value
cv = tinv(0.975, n1 - 1);

% Since the test statistic is barely smaller than the critical value, the
% null hypothesis is weakly accepted.

%% Two Sample t Test of the Logs of Data
% Since the log of the data appears much more Gaussian in the
% q-q plots, there is reason to believe that the two-sample t test 
% will preform better on them.

[h2, p2] = ttest2(log(clouds(:,1)),log(clouds(:,2)));

% The null hypothesis is now strongly rejected by the null hypothesis.

% The distribution of the data must be characterized before blindly
% applying a statistical test.