%% Example 8.15
% Use a permutation test on the Michelson speed of light data to assess
% whether the location parameters in each year are the same.

%% Importing Data
clear; close all; clf;
michelson1 = importdata('michelson1.dat');
michelson2 = importdata('michelson2.dat');

%% Permutation Test
n1 = length(michelson1);
n2 = length(michelson2);
n = n1 + n2;

michelson = [michelson1; michelson2];

s = sum(michelson1); % measure of the null hypothesis

sps = [];
b = 10000;
perm = zeros(b, n);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(n);
    if perm(i, :) == 1:n
        m = [m i];
    end
end

if ~isempty(m)
    perm(m',:) = []; % ensure original data are not included
end

perm = unique(perm, 'rows'); % ensure sampling without replacement

b = length(perm);

parfor i = 1:b
    michelsonp = michelson(perm(i,:)); % permuted version of the merged data
    sp = sum(michelsonp(1:n1));
    sps = [sps sp];
end

%% Permutation P-Value
pval = 2*min(sum(sps >= s) + 1, sum(sps < s) + 1)/(b + 1);

%% Plotting Permutation vs Gaussian
figure(1);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
x = 8.15*10^4:8.55*10^4;
pd = makedist('Normal','mu',mean(sps),'sigma',std(sps));
plot(x, pdf(pd,x),'r:','LineWidth',2)
xline(s,'k--')
hold off

ylabel 'Permutation Distribution'
legend 'Permutation Null Distribution' 'Gaussian Distribution with MLE' 'Test Statistic' ...
    'Location' 'southoutside'