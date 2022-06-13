%% Example 5.27
% The file geyser.dat contains 299 measurements of the waiting time in
% minutes between eruptions of Old Faithful in Yellowstone Nat. Park. A KDE
% shows that the distribution of data is bimodal. Fit a mixture Gaussian
% distribution to the data using the maximum likelihood estimator (MLE).

%% Calculating MLE parameters for mixture Gaussian distribution
clear; close all; clf
geyser = importdata('geyser.dat'); % geyser data

fun = @(lambda, data, cens, freq) -nansum(log(lambda(5)/(sqrt(2*pi)*lambda(3))...
    *exp(-((data - lambda(1))/lambda(3)).^2/2) + (1 - lambda(5))/(sqrt(2*pi)*...
    lambda(4))*exp(-((data - lambda(2))/lambda(4)).^2/2)));
% nansum() is used on purpose, to avoid overflow when the objective
% function is far from its maximum.

start = [55 80 5 10 .4]; % starting parameters
parm = mle(geyser, 'nloglf', fun, 'start', start, 'options', statset('MaxIter',300)); % mle parameters

%% Plotting KDE and mised Gaussian with MLE parameters
figure(1);

% Kernel Density Estimator (KDE)
[f,x] = ksdensity(geyser);
plot(x,f,'b')

% Mixture Gaussian Distribution
hold on
x = 20:0.1:140;
plot(x, parm(5)*normpdf(x,parm(1),parm(3)) + (1 - parm(5))*...
    normpdf(x,parm(2),parm(4)),'r')
hold off

ylabel 'Probability Density'
xlabel 'Waiting Time (min)'
legend('Kernel Density Estimator','Mixture Gaussian with MLE Paramters')
ylim([0 0.05])
title 'Old Faithful Eruption Waiting Times'

% The fitted distribution matches the means, but it is more sharply peaked,
% implying that the Gaussian distribution is not the appropriate model. We
% need something with flatter peaks...