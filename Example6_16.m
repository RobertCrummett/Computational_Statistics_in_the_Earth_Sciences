%% Example 6.16 - Barlett's M Test for Homogeneity of Variance
% This test was introduced by Barlett in 1937. It is more powerful than the
% F test, and applies to a population of independant variance estimates. It
% is also more robust to heteroskedasticity than the F test and is
% unbiased.

% I wrote a function to automatically computes the test statistic.

% For small samples, Barlett showed that his test statistic M is
% distributed according to chi squared (population size - 1)

% Because the null and alternate distirbutions cannot be fully
% characetrized, power calculations are not easy - however, the asymptotic
% power can be calculated from the noncentral chi squared distribution

%% Importing Data
clear; close all; clf;
data1 = importdata('michelson1.dat');
data2 = importdata('michelson2.dat');

%% Barlett's M Test
[M, sp2] = Bartlett(data1, data2);

pval = 2 * min(1 - chi2cdf(M,1), chi2cdf(M,1));

% This is weak evidence for the alternative hypothesis. This suggests that
% the data are nearly homoskedastic, unlike the F test.