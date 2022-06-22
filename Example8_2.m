%% Example 8.2
% For the Old Faithful data set in Examples 5.27 and 8.1, estimate the
% skewness along with its bias and standard error using the bootstrap.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');

%% Estimating Skewness with Bootstrapping Method
n = length(geyser);

rng default % this line can be removed to randomly sample multiple different ways

b = 1000; % number of bootstrap replicates
theta = skewness(geyser); % estimator

ind = unidrnd(n, n, b);
xboot = geyser(ind); % resample the data
thetab = skewness(xboot);

%% Results
meanb = mean(thetab); % bootstrap mean

baisb = mean(thetab) - theta; % bootstrap bias

stdb = std(thetab); % bootstrap standard error

bcb = 2*theta - mean(thetab); % bias corrected bootstrap estimate

%% Results from the MATLAB boostrp() Function
bmat = bootstrp(b, @skewness, geyser);

meanb2 = mean(bmat); % bootstrap mean

baisb2 = mean(bmat) - theta; % bootstrap bias

stdb2 = std(bmat); % bootstrap standard error
