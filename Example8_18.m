%% Example 8.18
% Use a permutation approach to test that the dispersion in the two
% Michelson speed of light data sets are the same after removing
% duplicates.

%% Importing Data
clear; close all; clf;
x1 = importdata('michelson1.dat');
x2 = importdata('michelson2.dat');

%% Removing Duplicates
x1 = unique(x1);
x1 = x1 - median(x1);
x2 = unique(x2);
x2 = x2 - median(x2);

%% Two Sample Permutation Test for Dispersion
x = [x1; x2];
n1 = length(x1);
n2 = length(x2);

s = sum(x1.^2)/n1 - sum(x2.^2)/n2; % test statistic

b = 10000;
perm = zeros(b, n1 + n2);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(n1 + n2);
    if perm(i,:) == 1:(n1 + n2)
        m = [m i];
    end
end

if ~isempty(m)
    perm(m',:) = []; % ensure original data are not included
end

perm = unique(perm,'rows'); % ensure sampling without replacement

b = length(perm);
sps = zeros(b, 1);

parfor i = 1:b
    xp = x(perm(i,:));
    sp = sum(xp(1:n1).^2)/n1 - sum(xp((n1+1):(n1 + n2)).^2)/n2;
    sps(i) = sp;
end

%% P-Value
pval = 2*min(sum(sps >= s) + 1, sum(sps < s) + 1)/(b + 1);

% The null hypothesis is accepted. The permutation distribution is slightly
% asymmetric.

%% Plotting KDE and Gaussian
figure(1);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
pd = makedist('Normal','mu',mean(sps),'sigma',std(sps));
xi2 = (-2*10^4):(2*10^4);
plot(xi2, pdf(pd, xi2), 'r:','LineWidth',2)
xline(s,'k--','LineWidth',2)
hold off

ylabel 'Permutation Distribution'
xlim([-2*10^4 2*10^4])

%% Repeating with Different Measure of Dispersion...
% An alternative measure of dispersion that is more robust is the
% difference of the sum of the absolute values of the deviations from the
% median.
x = [x1; x2];
n1 = length(x1);
n2 = length(x2);

s = sum(abs(x1))/n1 - sum(abs(x2))/n2; % test statistic

b = 10000;
perm = zeros(b, n1 + n2);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(n1 + n2);
    if perm(i,:) == 1:(n1 + n2)
        m = [m i];
    end
end

if ~isempty(m)
    perm(m',:) = []; % ensure original data are not included
end

perm = unique(perm,'rows'); % ensure sampling without replacement

b = length(perm);
sps = zeros(b, 1);

parfor i = 1:b
    xp = x(perm(i,:));
    sp = sum(abs(xp(1:n1)))/n1 - sum(abs(xp((n1+1):(n1 + n2))))/n2;
    sps(i) = sp;
end

%% P-Value
pval = 2*min(sum(sps >= s) + 1, sum(sps < s) + 1)/(b + 1);

%% Plotting KDE and Gaussian
figure(2);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
pd = makedist('Normal','mu',mean(sps),'sigma',std(sps));
xi2 = -80:0.1:80;
plot(xi2, pdf(pd, xi2), 'r:','LineWidth',2)
xline(s,'k--','LineWidth',2)
hold off

ylabel 'Permutation Distribution'
xlim([-80 80])

% My plots do not seem to match Dr. Chave's exactly... The problem is most
% likely on my end somewhere. But as I cannot presently see it, I will move
% on in confident expectation that the next Example will not be so trying.