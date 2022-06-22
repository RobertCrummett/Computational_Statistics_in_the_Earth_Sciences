%% Example 8.3
% Evaluate the effect of changing the number of bootstrap replicates for
% the skewness statistic using the geyser data.

%% Importing Data
clear; close all; clf;
geyser = importdata('geyser.dat');

%% Evaluate Bootstrap Replicates
rng default

b = [50 500 5000 50000];

figure(1);
for i = 1:length(b)
    bmat = bootstrp(b(i), @skewness, geyser);
    subplot(2,2,i)
    ksdensity(bmat)
    ylim([0 5])
    title(strcat(string(b(i))," Replicates"))
end
sgtitle 'KDEs of the bootstrap skewnesss for the Old Faithful data'

% The results indicate that a large number of samples results in a
% smoother, more symmetric distribution.s