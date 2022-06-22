%% Example 8.17
% Returning to the leukemia remission data from Example 7.15, use a
% permutation approach to assess whether a new drug improves survival.

%% Importing Data
clear; close all; clf;
remiss = importdata('remiss.dat');

%% Two-Sample Test for Paired Data
x = remiss(1:21, 2);
y = remiss(22:42, 2);

s = sum(x); % test statistic

data = [x; y];
n = length(data);
b = 100000;
perm = zeros(b,n);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(n);
    if perm(i,:) == 1:n
        m = [m i];
    end
end

if ~isempty(m)
    perm(m',:) = []; % ensure the original data are not included
end

perm = unique(perm, 'rows'); % ensure sampling without replacement

b = length(perm);
sps = zeros(b, 1);

parfor i = 1:b
    datap = data(perm(i,:));
    sp = sum(datap(1:n/2));
    sps(i) = sp;
end

%% Permutation P-Value
% Two sided, so multiply by 2
pval = 2*min(sum(sps >= s) + 1, sum(sps < s) + 1)/(b + 1);

%% Plotting KDE of Permutation Distribution
figure(1);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
pd = makedist('Normal','mu',mean(sps),'sigma',std(sps));
plot(150:400, pdf(pd, 150:400),'r:','LineWidth',2)
xline(s,'k--')
hold off

ylabel 'Permutation Distribution'
