%% Example 8.6
% Apply the bootstrap to get the double sided 95% confidence limits on
% skewness and kurtosis using bootstrap confidence intervals, t-strap
% confidence intervals, bootstrap percentile method, and bias-corrected
% and accelerated (BCa) method.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');
rng default % for consistent results

%% Bootstrap Confidence Interval on Skewness
bsskew = zeros(3, 2);

b = 10000; % number of replicates
bmat = bootstrp(b, @skewness, geyser);
xm = mean(bmat);
xs = std(bmat);

% Standard Confidence Interval
t = tinv(0.975, length(geyser) - 1);
bsskew(1,:) = [xm - t*xs, xm + t*xs];

% Bootstrap-t Intervals
zhat = (bmat - xm)./xs;
zhat = sort(zhat);
bsskew(2,:) = [xm + zhat(round(0.025*b))*xs, xm + zhat(round(0.975*b))*xs];

% Percentile Method
bmats = sort(bmat);
bsskew(3,:) = [bmats(round(b*0.025)), bmats(round(b*0.975))];

% BCa Method
ciskew = bootci(b, @skewness, geyser);

%% Bootstrap Confidence Interval on Kurtosis
bskurt = zeros(3, 2);

bmat = bootstrp(b, @kurtosis, geyser);
xm = mean(bmat);
xs = std(bmat);

% Standard Confidence Interval
t = tinv(0.975, length(geyser) - 1);
bskurt(1,:) = [xm - t*xs, xm + t*xs];

% Bootstrap-t Intervals
zhat = (bmat - xm)./xs;
zhat = sort(zhat);
bskurt(2,:) = [xm + zhat(round(0.025*b))*xs, xm + zhat(round(0.975*b))*xs];

% Percentile Method
bmats = sort(bmat);
bskurt(3,:) = [bmats(round(b*0.025)), bmats(round(b*0.975))];

% BCa Method
cikurt = bootci(b, @kurtosis, geyser);

% In general practice, the Percentile Method works the best when the numger
% of bootstrap replicates is large. It has better stability and convergance
% than either the Standard or T-Strap Methods. 

% For general applications, the BCa is the most consistent estimator.