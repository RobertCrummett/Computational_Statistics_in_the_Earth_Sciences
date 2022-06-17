%% Example 7.10
% Asses the fit of Gaussian distribution for the contaminated data of
% Example 7.8 and place confidence bounds on variance stabilized p-p plots
% and q-q plots of the data.

%% Data
clear; close all; clf;
rng default;
data = normrnd(0, 1, 1000, 1)';

% Contamination
data([1000, 999, 998, 997]) = [6.5 6.0 5.5 5];

%% P-P Plot
[h, p, ksstat, cv] = kstest(data);

% Uniform distribution
n = length(data);
u = ((1:n) - 0.5)/n;
r = 2/pi*asin(sqrt(u));

% P-P Plot
figure(1);
plot(r, 2/pi*asin(sqrt(normcdf(sort(data), mean(data), std(data)))), 'k+')
hold on
plot(r, 2/pi*asin(sqrt(min(1, sin(pi*r/2).^2 + cv))),'r--','LineWidth',1.5)
plot(r, 2/pi*asin(sqrt(max(0, sin(pi*r/2).^2 - cv))),'r--','LineWidth',1.5)
hold off

ylabel 'Transformed Data'
xlabel 'Sine Quantiles'
title 'Variance Stabilized P-P Plots with 95% Confidence Band from K-S Test'

%% Q-Q Plot
x = norminv(u, mean(data), std(data));

figure(2);
plot(x, sort((data - mean(data))/std(data)), 'k+')
hold on
plot(x, norminv(u + cv, mean(data), std(data)), 'r--','LineWidth',1.5)
plot(x, norminv(u - cv, mean(data), std(data)), 'r--','LineWidth',1.5)
hold off

ylabel 'Standardized Data'
xlabel 'Normal Quantiles'
title 'Q-Q Plot with 95% Confidence Band from K-S Test'

% In both plots, the widening of the confidence band does not allow for
% direct detection of the four upper outliers. However, in the Q-Q Plot,
% the slope has been tourqued so that is crosses the lower confidence
% bound. This would be reason to reject the null hypothesis that the data
% are Gaussian, or at the very least be a cause of suspicion.
