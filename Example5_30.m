%% Example 5.30
% Compute the Bonferroni confidence intervals for the parameters of the
% geyser data, example 5.27.

%% Import Data
clear; close all; clf;
data = importdata('geyser.dat');

%% MLE Parameters for misture Gaussian distribution

% Negative log likelihood function.
fun = @(lambda, data, cens, freq) - nansum(log(lambda(5)/(sqrt(2*pi)*lambda(3))*...
    exp(-((data - lambda(1))/lambda(3)).^2/2) + (1-lambda(5))/(sqrt(2*pi)*lambda(4))*...
    exp(-((data - lambda(2))/lambda(4)).^2/2)));

start = [55 80 5 10 .4]; % guess

% MLE parameters
parm = mle(data,['nloglf'],fun,'start',start,'Options',statset('MaxIter',300));

%% Asymptotic Covariance Matrix
acov = mlecov(parm, data, 'nloglf', fun);

xp = norminv(1 - 0.05/5);
sigma = sqrt(diag(acov)); % diagonal terms of acov are ignored, incorrectly
confBounds = [parm' - xp*sigma, parm' + xp*sigma];

% The confidence interval on the means is much tighter than the confidence
% intervals for the other intervals.