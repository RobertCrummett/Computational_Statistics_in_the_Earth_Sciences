%% Example 10.5
% Evaluate the equivalence of he covariance matrices for the global
% temperature data.

%% Importing Data
clear; clf; close all;
x = importdata('ZonAnnTs+dSST.dat'); 

% N and S hemisphere data in columns 3 and 4

%% Data
p = 4;
m = 13;
n = 10;
i = 1:n;

for j = 1:m
    xh(i, j, 1:p) = x((j - 1)*n + i, 8:(8 + p - 1));
end

sp = 0;
s1 = 0;

for i = 1:m
    % Summations
    s = cov(squeeze(xh(:, i, :)));
    sp = sp + (n - 1)*s; % pooled covariance matrix
    s1 = s1 + (n - 1)*log(det(s)); 
end

sp = sp/(m*(n - 1));
delta = (2*p^2 + 3*p - 1)/(6*(p + 1)*(m - 1))*(m/(n - 1) - 1/(m*(n - 1)));
u = ( 1- delta)*(m*(n - 1)*log(det(sp)) - s1);

pval = 2*min(chi2cdf(u, (m - 1)*p*(p + 1)/2), 1 - chi2cdf(u, (m - 1)*p*(p + 1)/2));

% The test accepts the null hypothesis that the covariance matrices for the
% 13 groups are the same. Consequently, it appears that the 10-year means
% of the data are different, but their covariance structures are similar,
% and hence MANOVA analysis in Example 10.2 is solid.