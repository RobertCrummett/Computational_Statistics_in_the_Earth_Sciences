%% Example 8.14
% Returning to Example 6.8, devise a permutation test for the null
% hypothesis that the number of days per year that an oceanographer spends
% at sea is 30 against an upper tail alternate. The MATLAB function
% randperm(n) produces a random permutation of the integers from 1 to n and
% will be used to obtain the permutations.

%% Data
clear; close all; clf;
x = [54 55 60 42 48 62 24 46 48 28 18 8 0 10 60 82 90 88 2 54];

%% Permutaion Test
% Test Statistic
x1 = sum(x - 30); 

n = length(x);

% Signs of the Differences
s = sign(x - 30);

b = 10000;
perm = zeros(b, n);
m = [];

parfor i = 1:b
    perm(i,:) = randperm(n);
    if perm(i,:) == 1:n
        m = [m i]
    end
end

if ~isempty(m)
    % Ensure original data are not included
    perm(m',:) = [];
end

% Ensuring Sampling Without Replacement
perm = unique(perm, 'rows');
b = length(perm);
sps = zeros([b 1]);

parfor i = 1:b
    % Permute the signs
    sp = s(perm(i, :));
    % Permutation Test Statistic
    x2 = sum(sp.*abs(x - 30));
    sps(i) = x2;
end

%% Exact P-Value
pval = (sum(sps >= x1) + 1)/(b + 1);

% The null hypothesis is weakly rejected.

%% Plotting
figure(1);
[f, xi] = ksdensity(sps);
plot(xi, f, 'b-','LineWidth',2)
hold on
pd = makedist("Normal","mu",mean(sps),"sigma",std(sps));
plot(-100:400, pdf(pd, -100:400),'r:','LineWidth',2)
xline(x1,'g--','LineWidth',2)
ylim([0 8*10^(-3)])
hold off

ylabel 'Permutation Distribution';
legend('KDE of Permutation Distirbution','Gaussian PDF with MLE params',...
    'Test Statistic','Location','nw')

% The region of rejection is to the right of the test statistic in the
% figure.