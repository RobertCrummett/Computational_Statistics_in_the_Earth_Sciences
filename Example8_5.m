%% Example 8.5
% Applt the smoothed bootstrap to the interquartile range of the geyser
% data.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');

%% Smoother Bootstrapping Method
rng default

h = [0.707, 0.866, 1, 1.414];
b = 50000; % numebr of replicates
bmat = bootstrp(b, @iqr, geyser);

figure(1)
for i = 1:length(h)
    x = bmat + normrnd(bmat, h(i).^2*ones(size(bmat)));
    subplot(2, 2, i)
    ksdensity(x)
    title(strcat("Smoothing Bandwidth = ",string(h(i))))
    ylim([0 0.15])
end
sgtitle(["KDE of 50000 bootstrap replicates for the IQR";...
    "Smoothed by Gaussian kernel with differing bandwidths"])

% Smoothing bandwidth of about 1 seems appropriate.