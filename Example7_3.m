%% Example 7.3
% The data set contains 62 measurements of the time in days between
% successive earthquakes taken from Hand et al. (1994). Assess the
% goodness-of-fit of the data to the exponential distribution.

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% Plot Data verse Exponential Distribution
figure(1);
[f, x] = ecdf(quakes);
plot(x, f, 'b')
hold on
plot(x, expcdf(x, mean(quakes)),'r')
hold off

legend('Empirical CDF','Exponential CDF','Location','se')
xlabel 'Time Between Earthquakes (Days)';
ylabel 'Cumulative Probability';
title 'The Empirical and Exponential CDFs for the Earthquake Interval Time Data'

% The exponential cdf has the maximum likelihood estimate paramter.

%% Histogram
figure(2);
h = histogram(quakes, 'NumBins', 7, 'Normalization', 'count');
title 'Histogram of Earthquake Count Data'
xlabel 'Time Between Earthquakes (Days)';
ylabel 'Count';

%% Likelihood Ratio Testing
% Observed
o = h.Values;
p = zeros([1 length(o)]);

for i = 1:length(o)
    p(i) = integral(@(x) exppdf(x, mean(quakes)), h.BinEdges(i), h.BinEdges(i) + h.BinWidth);
end

% Expected
e = p*length(quakes);

% Test Statistic
loglambda = 2*sum(o.*log(o./e));

% Critical Value
cv = chi2inv(0.975, length(o) - 2);

% P-Value
pval = 1 - chi2cdf(loglambda, length(o) - 2);

% The null hypothesis that the data is exponential is accepted. However,
% the number of data per bin is far from even. They will be adjusted to
% make that more true at the cost of fewer bins.

%% Smaller Bins Histogram
figure(3);
edges = [ 0 75 200 400 600 1000 1960];
h = histogram(quakes, edges);
title 'Histogram of Earthquake Count Data'
xlabel 'Time Between Earthquakes (Days)';
ylabel 'Count';

%% Likelihood Ratio Testing with Fewer, More Even, Bins
% Observed
o = h.Values;
p = zeros([1 length(o)]);

for i = 1:length(o)
    p(i) = integral(@(x) exppdf(x, mean(quakes)), h.BinEdges(i), h.BinEdges(i + 1));
end

% Expected
e = p*length(quakes);

% Test Statistic
loglambda2 = 2*sum(o.*log(o./e));

% Critical Value
cv2 = chi2inv(0.975, length(o) - 2);

% P-Value
pval2 = 1 - chi2cdf(loglambda, length(o) - 2);

% The previous conclusion remains unchanged.