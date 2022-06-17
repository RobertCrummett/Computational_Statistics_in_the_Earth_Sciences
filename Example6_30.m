%% Example 6.30
clear; close all; clf;
% Suppose 10 independant hypothesis tests are computed, leading to the
% following set of ordered p-values:
pval = [0.00013 0.00333 0.00579 0.00956 0.01345 0.24689 0.35487 0.55698 0.65345 0.99111];

% Apply a Bonferroni and Benjamini-Hochberg multiple hypothesis test at the
% 0.05 level.

%% Bonferroni
A = 0.05;
N = length(pval);

% Threshold
alpha = A/N;
Rej = nnz(find(pval <= alpha));

% Thus, the first two p-values are rejected.

%% Benjamini-Hochberg
[RejTh, Index] = Benjamini(pval);

% Thus, the first five p-values are rejected.

