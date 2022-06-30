%% Example 10.11
% Use empirical orthogonal function regression to evaluate the predictors
% of multicollinearity by applying a rank-reduced solution for stability
% and then by applying Wilk's lambda test for multivariate regression to
% test the result.

%% Importing Data
clear; clf; close all;
sake = importdata('sake.dat');

%% Data
y = sake(:, 1:2);
y = y - repmat(mean(y), length(y), 1);
x = sake(:, 3:10);
x = x - repmat(mean(x), length(x), 1);

[n, m] = size(x);

%% Rank-Reduced Solutions
[u, s, v] = svd(x);

for i = m:-1:2
    beta = v(:,1:i)*inv(s(:,1:i)'*s(:,1:i))*s(:,1:i)'*u'*y;
    e = (y - u*s(:,1:i)*v(:,1:i)'*beta)'*y;
    h = beta'*v(:,1:i)*s(:,1:i)'*u'*y;
    p(i) = WilksLambda(h, e, 2, i, n - i - 1);
end

% The p-values are all well above 0.05, and hence there is no need for
% regresion even when only two eigenvalues are retained. The problem with
% these data is not multicollinearity. Further, the result is not changed
% if the variables x amd y are standardized so that their principal
% components are those from the correlation matrix.