%% Example 7.11
% For the earthquake data, compute the Anderson-Darling (A-D) test for the
% null hypothesis that the data are exponential against the alternate that
% they are not.

% The Anderson-Darling test is more powerful than the Kolmogorov-Smirnov
% test. It is sensitive to the distribution tailsm but has the disadvantage
% that the critical values depend on the target distribution.

%% Importing Data
clear; close all; clf;
quakes = importdata('quakes.dat');

%% A-D test
quakes = sort(quakes/mean(quakes));
N = length(quakes);
i = 1:N;
u1 = expcdf(quakes);
u2 = 1 - expcdf(quakes(N - i + 1));
a2 = - N - sum((2*i - 1)'.*(log(u1) + log(u2)))/N;

% The critical value at alpha = 0.05 is 1.326, so the null hypothesis is
% accepted.

%% MATLAB adtest() Function
[h, p, adstat, cv] = adtest(quakes,'Distribution','exp');

% The p-value is high, so the null hypothesis is strongly supported. So
% the data are exponential.