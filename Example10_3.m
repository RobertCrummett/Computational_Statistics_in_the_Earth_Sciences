%% Example 10.3
% Determine whether there is any difference in the skull sizes over time.

%% Importing Data
clear; clf; close all;
skulls = importdata('skulls.dat');

%% Data
p = 4;
m = 5;
n = 30;

% Preallocation
x = zeros(n, m, p);

for i = 1:n
    for j = 1:m
        x(i, j, 1:p) = skulls((j - 1)*n + i, 1:p);
    end
end

% Preallocation
e = zeros(p,p);
h = zeros(p,p);

for i = 1:p
    for j = 1:p
        h(i, j) = sum(sum(x(:, :, i)*(eye(m) - ones(m,m)/m)*x(:, :, j)'))/n;
        e(i, j) = trace(x(:, :, i)*x(:, :, j)') - sum(sum(x(:, :, i)*x(:, :, j)'))/n;
    end
end

% Wilk's Lambda Test
[pval, lambda] = WilksLambda(h, e, p, m - 1, n*(m - 1));

% The p-value is very small. Therefore the null hypothesis that there is no
% difference in skull size is rejected.

%% ANOVA Tables
% ANOVA tables for each variable over time will be examined to assess the
% univariate relations among the data.

for j = 1:p
    [pp, anovatab] = anova1(x(:, :, j));
end

% Only the last variable (nasal height of skull) accepts the univariate
% null hypothesis, althought the second one only weakly rejects it.
% If Bonferroni approach is employed, the second and last variables accept
% the null hypothesis.
% Same with Benjamini-Hochberg method implementation.
% The remaining variables all demonstrate significant change through time.