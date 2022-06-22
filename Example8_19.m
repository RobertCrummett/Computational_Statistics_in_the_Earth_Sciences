%% Example 8.19
% Use the jackknife to estimate the skewness and its standard error for the
% Old Faithful data set.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');

%% Jackknife
n = length(geyser);

skew = skewness(geyser);

skewt = zeros(1, n);
for i = 1:n
    geyt = geyser;
    geyt(i) = []; % keave i-th point out
    skewt(i) = skewness(geyt);
end

jkskew = mean(skewt); % exact same value as using the skewness() function on original data
jkskewvar = sqrt((n -1)/n*sum((skewt - mean(skewt)).^2));

%% MATLAB jackknife() Function
jackstat = jackknife(@skewness, geyser);

jkskew2 = mean(jackstat);
jkskewvar2 = sqrt((n -1)/n*sum((skewt - mean(skewt)).^2));

%% Example 8.20
% Now estimate the median and its standard error from the geyser data using
% the jackknife and bootstrap methods.

%% Bootstrapping Median of Geyser Data
b = 100000; % repititions

idx = unidrnd(n, n, b);
boot = geyser(idx);
bmat = median(boot);

% Standard Error for Bootstrap
bmed_se = sqrt((n - 1)/n*sum(bmat - mean(bmat)).^2);

%% Jackknife Median of Geyser Data
medt = zeros(1, n);

for i = 1:n
    geys = geyser;
    geys(i) = [];
    medt(i) = median(geys);
end

% Standard Error for Jackknife
jkmed_se = sqrt((n - 1)/n*sum(medt - mean(medt)).^2);

% Since the jackknife method breaks down estimators based on order
% statistics, the mdeian fails, because it is not a smooth statistic.

