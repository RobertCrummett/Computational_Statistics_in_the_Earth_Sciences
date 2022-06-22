%% Example 8.9
% Evaluate the null hypothesis that the days an oceanographer spends at sea
% per year is 30 days against the alternative that it is not using 99, 999,
% 9999, and 99999 bootstrap replicates.

%% Data
clear; close all; clf;
data = [54 55 60 52 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];
rng default

%% Bootstrap Hypothesis Testing
b = [99 999 9999 99999];
that = (mean(data) - 30)./std(data, 1); % Sample test statistic against the hypothesized value
pval = zeros([length(b) 1]); % preallocation

for i=1:length(b)
    bmat = bootstrp(b(i), @(t) (mean(t) - mean(data))./std(t,1), data);
    % Two-Tailed P-Value
    pval(i) = 2*min(sum(bmat >= that), sum(bmat < that))./b(i);
end

% The test rejects the null hypothesis that the oceanographer spends 30
% days at sea, although not strongly.

%% Plotting KDE of Test Statistic
figure(1);
ksdensity(bmat)
xline(that,'r--')
ylim([0 2])
title 'KDE for 99,999 bootstrap replicates testing H_0: \mu* = 30'
dc = [double('t'),770];
legend('KDE',char(dc),'Location','nw')
