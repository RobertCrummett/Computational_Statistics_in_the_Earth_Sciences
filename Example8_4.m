%% Example 8.4
% Apply the bootstrapping method to the interquartile range and explain the
% result.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');

%% Bootstrapping Method for IQR
rng default

b = [50 500 5000 50000];

figure(1);
for i = 1:length(b)
    bmat = bootstrp(b(i), @iqr, geyser);
    subplot(2,2,i)
    ksdensity(bmat)
    title(strcat(string(b(i))," Replicates"))
end
sgtitle 'Bootstrapping Method for Interquartile Range'

% As the number of replicates increases, multiple peaks in the distribution
% are resolved. This occurs because order statistics are concentrated on
% the sample values that are discrete unless the data set is very large.

% Thus, as the number of replicates rises, the bootstrapping method becomes
% more discrete, and the naive nonparametric bootstrapping approach fails.