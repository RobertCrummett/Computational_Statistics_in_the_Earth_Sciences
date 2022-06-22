%% Example 8.11
% Evaluate the null hypothesis that the data are exponential and Rayleigh
% by testing against random draws with MLE parameters.

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% Test Whether Two Samples Are Drawn From the Same Distribution
rng default
rnd = exprnd(mean(quakes), size(quakes));
data = [quakes; rnd]';

b = [99 999 9999 99999];

that = mean(quakes) - mean(rnd);

%% Bootstrapping, Testing Difference of Means of Exponential Distribution
% Preallocation
rng default
pval = zeros([1 length(b)]);

for i = 1:length(b)
    % Bootstapping Random Samples
    bmat = bootstrp(b(i), @(x) (mean(x(1:length(quakes))) - mean(x(length(quakes) + 1:2*length(quakes)))), data);
    % P-Values
    pval(i) = 2*min(sum(bmat > that), sum(bmat <= that))./b(i);
end

% The test statisitc is apparently 121.873, but I have no idea how Dr.
% Chave reached that conclusion. that is plotted in the figure.

%% Plotting Bootstrap Data for Difference of Means
figure(1);
ksdensity(bmat)
xline(that,'r--')
title 'KDE of the Difference of Means Bootstrap Data'
legend('KDE','test statistic','Location','ne')

%% Bootstrapping, Testing Difference of Medians of Exponential Distribution
% Test Statistic
rng default
that_med = median(quakes) - median(rnd);

% Preallocation
pval_med = zeros([1 length(b)]);

for i = 1:length(b)
    % Bootstapping Random Samples
    bmat = bootstrp(b(i), @(x) (median(x(1:length(quakes))) - median(x(length(quakes) + 1:2*length(quakes)))), data);
    % P-Values
    pval_med(i) = 2*min(sum(bmat > that_med), sum(bmat <= that_med))./b(i);
end

% The null hypothesis that the data are exponential is accpeted once again.

%% Plotting Bootstrap Data for Difference of Medians
figure(2);
ksdensity(bmat)
xline(that_med,'r--')
title 'KDE of the Difference of Medians Bootstrap Data'
legend('KDE','test statistic','Location','nw')

%% Bootstrapping, Testing Difference of Means of Rayleigh Distribution
% MLE parameter
mlerayl = sqrt(1/length(quakes)*(sum(quakes.^2/2)));

% Random Draws
rng default
rnd = raylrnd(mlerayl, size(quakes));
data = [quakes; rnd]';

that_rayl = mean(quakes) - mean(rnd);

% Preallocation
pval_rayl = zeros([1 length(b)]);

for i = 1:length(b)
    % Bootstapping Random Samples
    bmat = bootstrp(b(i), @(x) (mean(x(1:length(quakes))) - mean(x(length(quakes) + 1:2*length(quakes)))), data);
    % P-Values
    pval_rayl(i) = 2*min(sum(bmat > that_rayl), sum(bmat <= that_rayl))./b(i);
end

% The null hypothesis that the data are Rayleigh is strongly rejected. The
% bootstrap data is long tailed.

%% Plotting Bootstrap Data for Difference of Medians
figure(3);
ksdensity(bmat)
xline(that_rayl,'r--')
title 'KDE of the Difference of Mean Bootstrap Data for Rayleigh Distribtuion'
legend('KDE','test statistic','Location','ne')