%% Example 6.19
% Determine the number of data with the observed mean and standard
% deviation that would be required to achieve a power of 0.90 at a
% significance level of 0.05 for a two sided test with population means of
% 30 and 40.

%% Data
clear; close all; clf;
x = [54 55 60 42 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];

%% Sample Size Necessary for Different Population Means
% Population means tested are 30 and 40.
n30 = sampsizepwr('t', [mean(x) std(x)], 30);
n40 = sampsizepwr('t', [mean(x) std(x)], 40);

% The smaller the distance between the population mean and the sample mean,
% the larger the sample size needed to resolve the difference.

%% Power of t Test with Sample Size 20, 30, and 100
% Population mean is 30.
b20 = sampsizepwr('t', [mean(x) std(x)], 30, [], 20);
b30 = sampsizepwr('t', [mean(x) std(x)], 30, [], 30);
b100 = sampsizepwr('t', [mean(x) std(x)], 30, [], 100);

% Clearly, asking 30 oceanographers rather than 20 would have been a much
% more powerful test.