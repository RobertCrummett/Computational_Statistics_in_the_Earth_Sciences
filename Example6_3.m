%% Example 6.3 - 6.5
% Type two errors depend on the true mean - therefore, errors must be
% evaluated over a range of values.

%% Calculating Beta Values for Different Means
clear; close all; clf;

mualt = 10:0.1:30;
cv = norminv(0.95); % one-sided critical value
sigma = 0.6;
ct = cv + 20; % critical value + mean
ctv = ct * ones(size(mualt));

beta = normcdf(ctv, mualt, sigma); % probabiltity of type 2 error
power = 1 - beta; % power
p = 0.05 * ones(size(mualt));

%% Plotting the Power Curve
figure(1);
plot(mualt, power, mualt, p)
xlabel 'True Mean'; 
ylabel 'Power'; 
title 'Power Curve for the Climate Example';
legend('Power curve','\alpha = 0.05','Location','nw')

%% P-Value Calculated
mu = 20;
sig = 0.6;
xbar = 22;
zobs = (xbar - mu)/sig;
pval = 1  - normcdf(zobs, 0, 1); 

% The very small p-value suggests that the alternative hypothesis is
% correct. Combined with the reasonably high power of the test, (about 72%
% for the observed mean of 22) the null hypothesis that the observed mean
% is reasonable can be rejected confidently.