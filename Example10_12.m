%% Example 10.12
% There is one response variable, six predictor variables, and 16
% observations in the seven collumns of the data file. The predictors are
% highly incommensurable, and hence the data will be standardized for
% analysis.

%% Importing Data
clear; clf; close all;
data = importdata('longley.dat');

%% Data Standardization
y = data(:, 1);
y = y - repmat(mean(y), length(y), 1);
y = y./repmat(std(y), length(y), 1);
x = data(:, 2:7);
x = x - repmat(mean(x), length(x), 1);
x = x./repmat(std(x), length(x), 1);

[n, p] = size(x);
w = ones(size(y));

[b, bse, t, pval] = Treg(y, x, w);

% For the first and fifith predictors, the standard errors (bse) are
% substantially larger than the coefficents, and the corresponding
% t-statistics are very small. A test that the coefficents are zero rejects
% only the third, fourth, and sixth predictors. The problem is that the
% predictors are highly correlated.

disp(corr(x)) % wow!

%% Singular-Value Decomposition
[u, s, v] = svd(x);

percvar = diag(s.^2)/sum(diag(s.^2));

% The first three eigenvalues account for 99.7% of the variance and will be
% retained for further analysis.

%% EOF Regression
bp = v(:,1:3)*inv(s(:,1:3)'*s(:,1:3))*s(:,1:3)'*u'*y;
yhat = u*s(:,1:3)*v(:,1:3)'*bp;
sigma2 = (y - yhat)'*(y - yhat)/(n - p);
bpse = sqrt(sigma2*diag(v(:,1:3)*inv(s(:,1:3)'*s(:,1:3))*v(:,1:3)'));
tp = bp./bpse;
pvalp = 2*(1 - tcdf(abs(tp), n - p));

% The p-values indicate that the regression coefficents are required for
% all but the fourth predictor. And the acceptance of the null hypothesis
% for that predictor is weak. The regression coefficents have chagned
% dramatically, indicating the influence of multicollinearity in the data.