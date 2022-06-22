%% Example 8.8
% From the earthquake data in Example 7.3, compute the skewness and
% kurtosis along with bootstrap precentile and BCa confidence intervals.

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');
rng default % reproduceability

%% Skewness and Kurtosis
skew = skewness(quakes);
kurt = kurtosis(quakes);

%% Bootstrap Confidence Interval Estimation of Skewness and Kurtosis
b = 10000;

% Skewness Confidence Interval Estimation
% Percentile
bootci(b, {@skewness, quakes},'type','per')
% BCa
bootci(b, @skewness, quakes)

% Kurtosis Confidence Interval Estimation
% Percentile
bootci(b, {@kurtosis, quakes},'type','per')
% BCa
bootci(b, @kurtosis, quakes)

% The differences between the estimators is substantial! The difference is
% not due to the bias correction of BCa. The difference is accounted for by
% recognizing that the distributions of the skewness and kurtosis are
% themselves skewed.

% In this circumstance, the percentile confidence interval will be biased,
% whereas the BCa interval will offer better preformance.

%% Plotting KDE of Bootstrap Distributions

figure(1);
subplot(1,2,1)
bmat = bootstrp(b, @skewness, quakes);
ksdensity(bmat);
title 'Skewness PDF'
axis square
ylim([0 1.4])

subplot(1,2,2)
bmat = bootstrp(b, @kurtosis, quakes);
ksdensity(bmat)
title 'Kurtosis PDF'
axis square
ylim([0 0.35])

sgtitle 'KDE for PDF for the bootstrap distributions with 10000 replicates'
