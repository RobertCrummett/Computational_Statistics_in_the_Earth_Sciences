%% Example 8.16
% Apply the permutation test to evaluate the null hypothesis that seeded
% and unseeded clouds produce the same amount of rain using the data from
% Example 6.12

%% Importing Data
clear; close all; clf;
data = importdata('cloudseed.dat');

%% Permutation Test
n = length(data);
cloud = [data(:,1); data(:,2)];
s = mean(data(:,1)) - mean(data(:,2));

b = 10000;
perm = zeros(b, 2*n);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(2*n);
    if perm(i,:) == 1:2*n
        m = [m i];
    end
end

if ~isempty(m)
    perm(m',:) = []; % ensure original data are not included
end

perm = unique(perm, "rows"); % ensure sampling without replacement
b = length(perm);
sps = zeros([b 1]);

parfor i = 1:b
    cloudp = cloud(perm(i,:));
    sp = mean(cloudp(1:n)) - mean(cloudp(n+1:2*n));
    sps(i) = sp;
end

%% Permutation P-Value
pval = 2*min(sum(sps >= s) + 1, sum(sps < s) + 1)/(b + 1);

% The null hypothesis that the seeded and unseeded clouds produce the same
% amount of rain is rejected, although only weakly.

%% Plotting Permutation Distribution
figure(1);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
pd = makedist('Normal','mu',mean(sps),'sigma',std(sps));
xi2 = -500:500;
plot(xi2, pdf(pd, xi2),'r:','LineWidth',2)
xline(s,'k--')
hold off

ylabel 'Permutation Distribution'
ylim([0 3*10^(-3)])
xlim([-500 500])
title 'Permutation Distribution vs Gaussian Distribution with MLE       '
