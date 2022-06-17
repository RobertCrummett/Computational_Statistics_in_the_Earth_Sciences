%% Example 6.25
% The mean monthly wind direction for a high-latitude region using the
% National Center for Environmental Prediction product is contained in the
% file wind.dat for 2004 - 2005. The monthly means will average out
% periodic phenomena. The posulated mean over the 2 year interval is 240
% degrees. Derive the likelihood ratio test for the null hypothesis nu =
% nu* and kappa > 0, verses the alternative nu does not equal nu* and kappa
% > 0 for the Von Mises Distribution.

%% Import Data
clear; close all; clf;
wind = importdata('wind.dat');
figure(1);
polarhistogram(wind,'NumBins',25)
title 'Rose diagram for the wind direction data'

%% Likelihood Ratio Testing
% For the Von Mises distribution.
s = sum(sind(wind));
c = sum(cosd(wind));
r = sqrt(s^2 + c^2);

% Maximum Likelihood Estimates for nu
nuhat = atan2d(s, c);
nustar = 240;

% Maximum Likelihood Estimates for kappa
N = length(wind);
khat = fzero(@(x) besseli(1, x)/besseli(0, x) - (r/N), 1);
khathat = fzero(@(x) besseli(1, x)/besseli(0, x) - ((c*cosd(nustar) + s*sind(nustar))/N), 1);

% Test Statistic
test =  2*khat*sum(cosd(wind - nuhat)) + 2*N*log(besseli(0, khathat)) -...
    2*khathat*sum(cosd(wind - nustar)) - 2*N*log(besseli(0, khat));

% Critical Value
cv = chi2inv(0.975, 1);

% P-Value
pval = 1 - chi2cdf(test, 1);

% The data is not rejecting the null hypothesis when nustar = 240.

%% Plotting Power Curves
%Variance Estimate
sigma2 = 1./(r*khat);

clear nustar
nustar = 230:0.1:250;
delta = sind(nuhat - nustar).^2/sigma2;
xhi = chi2inv(0.975, 1);
xlo = chi2inv(0.025, 1);
power = 1 - ncx2cdf(xhi, 1, delta) + ncx2cdf(xlo, 1, delta);

figure(2);
plot(nustar, power,'m','LineWidth',1.5)
xlabel 'Hypothesized Mean'; 
ylabel 'Power';
ylim([0.05 0.25])
title 'Power of the Wind Direction Likelihood Ratio Test for Von Mises Distribution';

% The test has low power, and hence a low probability of rejecting the null
% hypothesis when it is false.