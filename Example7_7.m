%% Example 7.7, 7.9
% Estimate the fit of the earthquake data of Example 7.3 to a Weibull
% distribution using a two sided Kolmogorov-Smirnov test.

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% MLE Parameters for Weibull CDF
parmhat = wblfit(quakes);

% Sorting Data
quakes = sort(quakes);

cdf = [quakes'; wblcdf(quakes', parmhat(1), parmhat(2))]';

%% K-S Test
[h, p, ksstat, cv] = kstest(quakes, cdf);

% The null hypothesis is accepted at the 0.05 level, and the p-value is
% large.

%% Plotting Weibull CDF against  Empirical CDF
figure(1);
[f, x] = ecdf(quakes);
plot(x, f, 'b')
hold on
plot(x, wblcdf(x, parmhat(1), parmhat(2)),'r')
hold off

legend('Empirical CDF','Weibull CDF with MLE parameters','Location','se')
xlabel 'Time Between Earthquakes (Days)';
ylabel 'Cumulative Probability';
title 'The Empirical and Weibull CDFs for the Earthquake Interval Time Data'

%% Plotting Confidence Bounds on the Weibull CDF
cv = 1.3581; % c_{0.05}, approximate critical value

figure(2);
plot(x, wblcdf(x, parmhat(1), parmhat(2)), 'k')
hold on
plot(x, max(0, f - cv/sqrt(length(quakes))),'r--')
plot(x, min(1, f + cv/sqrt(length(quakes))),'r--')
hold off

legend('Weibull CDF with MLE parameters','C_{0.05} Confidence Bounds','Location','se')
xlabel 'Time Between Earthquakes (Days)';
ylabel 'Cumulative Probability';
title 'Weibull CDF with Confidence Bounds for the Earthquake Interval Time Data'
