%% Example 8.13
% Remove the bias from the K-S test applied to contaminated Gaussian data
% in Example 7.8

%% Data
clear; close all; clf;
rng default
data = normrnd(0, 1, 1000, 1);

% Contamination
data([1000 999 998 997]) = [6.5 6 5.5 5];

%% K-S Test
[h, p, ksstat] = kstest(data);

par1 = mean(data);
par2 = std(data, 1);

b = 99999;
sps = zeros([1 length(b)]);

parfor i = 1:b
    draw = sort(normrnd(par1, par2, 1000, 1));
    cdf = [draw normcdf(draw, par1, par2)];
    [~,~,ks1] = kstest(draw, cdf);
    sps(i) = ks1;
end

pval = 2*min(sum(sps >= ksstat) + 1, sum(sps < ksstat) + 1)/(b + 1);

% The downward bias of the K-S statistic is evident. The conclusion remains
% unchanged however.