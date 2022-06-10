%% Example 4.10
% Obtain 10000 random draws from the arcsine distribution, and compute the
% pdf with support that is unbounded and confined to [-1 1]. The arcsine
% distribution si the beta distribution with parameters 1/2 and 1/2.
% Repeat this for the cdf.

%% Random draws
clear; close all; clf;
x = betarnd(0.5, 0.5, 10000, 1);

%% KDE
[f, xi, bw] = ksdensity(x);

%% Plot KDE and beta distribution pdf
% Plot KDE without support first
figure(1);
plot(xi,f,'b--')
hold on

% Plot KDE with support from [0 1]
[f, xi] = ksdensity(x, 'support', [0, 1]);
plot(xi, f, 'r')

% Plot beta(x,1/2,1/2), the arcsine distribution
plot(xi, betapdf(xi, 0.5, 0.5), 'k')
hold off
ylim([0 5]) %necessary to see what is going on
legend('no support','support from [0 1]','arcsine distribution')
title 'Random draws from arcsine distribution fit by a kernel density estimator'

%% Plot KDE cdf and beta distribution cdf
% Plot KDE without support first
[f, xi] = ksdensity(x,'function','cdf');
figure(2);
plot(xi, f, 'b--')

% Plot KDE with suppport from [0 1]
[f, xi] = ksdensity(x, 'support', [0 1], 'function', 'cdf');
hold on
plot(xi, f, 'r-')

% Plot betacdf(x, 1/2, 1/2), the cdf of the arcsine distribution
plot(xi, betacdf(xi, 0.5, 0.5), 'k-')
hold off
legend('no support','support from [0 1]','arcsine cdf','Location','nw')
title 'Random draws from arcsine cdf fit by a kernel density estimator'
