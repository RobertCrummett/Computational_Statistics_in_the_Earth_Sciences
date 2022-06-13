%% Example: Distribution of the Sample Median
% In this case, for a normal parent distribution.

%% Expected Value of Sample Median
clear; close all; clf;
N = 101; % sample size

fun = @(x, N) x.*normcdf(x).^floor(N/2).*(1 - normcdf(x)).^(N - floor(N/2) - 1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2));
muhat = integral(@(x) fun(x,N), -Inf, Inf); % expected value

% For different sample sizes, the expected value is floating point zero.
% -- For odd values of N only! --
% This suggests that the sample median is an unbaised estimator.

%% Variance of Sample Median
N = 1001; % sample size

fun = @(x, N) x.^2.*normcdf(x).^floor(N/2).*(1 - normcdf(x)).^(N - floor(N/2) - 1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2));
sig2hat = integral(@(x) fun(x,N), -Inf, Inf);

% Presuming that the variance goes like alpha*sigma^2/N, where alpha is a
% constant, the estimates of alpha converge to pi/2
% So the variance of the sample median is 60% larger than the sample mean
% for a normal parent

%% Plotting PDF of Sample Median from Gaussian parent distribution
x = -2:0.01:2; % range

figure(1);
N = 3;
plot(x, normcdf(x).^floor(N/2).*(1 - normcdf(x)).^(N - floor(N/2) - 1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2)),'b')
hold on
N = 11;
plot(x, normcdf(x).^floor(N/2).*(1 - normcdf(x)).^(N - floor(N/2) - 1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2)),'r')
N = 101;
plot(x,normcdf(x).^floor(N/2).*(1 - normcdf(x)).^(N - floor(N/2) - 1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2)),'g')
hold off

ylabel 'Median PDF'
ylim([0 4])
title 'Sampling distribution for the median for a standard normal parent distribution'
legend('N = 3','N = 11','N = 101','Location','ne')
