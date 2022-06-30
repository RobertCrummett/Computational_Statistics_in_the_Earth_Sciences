%% Example 10.4
% Apply the sphericity test to determine if the variables are correlated.

%% Importing Data
clear; clf; close all;
skulls = importdata('skulls.dat');

%% Sphericity Test
p = 4;
m = 5;
n = 30;

s = cov(skulls);
u = p^p*det(s)/trace(s)^p;
test = -n*m*log(u);

pval = 2*min(chi2cdf(test, p*(p + 1)/2 - 1), 1 - chi2cdf(test, p*(p + 1)/2 - 1));

% Since the p-value is zero, the sphericity test rejects the null hypothess
% that the data are uncorrelated. Thus, a MANOVA approach has to be used.