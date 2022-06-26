%% Example 9.1
% First, a regression analysis will be preformed on the data, including a
% constant term, and the results will be evaluated using ANOVA, t, runs,
% and Durbin-Watson tests.

%% Importing Data
clear; close all; clf;
hubble = importdata('hubble.dat');

%% Data Preparation
y = hubble(:, 2);
x = [hubble(:,1) ones(size(y))];
w = ones(size(y));

%% ANOVA
a = Anovar(y, x, w);

%% T-Table for Linear Regression
[b, bse, t, tprob] = Treg(y, x, w);

%% Run Test
% Residuals
r = y - x*b;
r2 = a(1,2)/a(3,2);

[h, p] = runstest(r);

%% Durbin-Watson Test
pd = dwtest(r, x);
% Not the same result as Dr. Chave!

% The D-W tests indicate that the regression residuals are random and
% uncorrelated.

%% Correlation Coefficent R2
Q2 = a(2,2)/a(3,2);

% Correlation Coefficent
R2 = 1 - Q2;

%% Intercept May Be Neglected for Further Analysis...
% Since the p-value for the intercepts accept the null hypothesis that they
% are zero, we can neglect the column of ones that give us intercepts,
% although it is strongly recommended we use it.

x = hubble(:,1);

%% ANOVA
a = Anovar(y, x, w);

%% T-Table for Linear Regression
[b2, bse, t, tprob] = Treg(y, x, w);

%% Run Test
r = y - x*b2;
r2 = a(1,2)/a(3,2);

[h, p] = runstest(r);

%% Durbin-Watson Test
pd = dwtest(r, x);

% Both the Durbin-Watson test and the runs test continue to show that the
% residuals are random and uncorrelated.

%% Correlation Coefficent R2
Q2 = a(2,2)/a(3,2);

% Correlation Coefficent
R2 = 1 - Q2;

% Slightly smaller than before.

%% Plotting Hubble Data with Linear Regression
figure(1);
plot(hubble(:,1), hubble(:,2), 'k+', 'LineWidth', 1.5)
hold on
plot(sort(hubble(:,1)), b(2) + b(1)*sort(hubble(:,1)),'r--', 'LineWidth', 2)
plot(sort(hubble(:,1)), b2*sort(hubble(:,1)),'b--', 'LineWidth', 2)
hold off

legend("Hubble's Data", "Lin. Reg. with Intercept","Lin. Reg. without Intercept",'Location','se')
title 'The Hubble Data with Best-Fit Linear Regression'
xlabel 'Distance (Mp)'
ylabel 'Recession Velocity (km/s)'
ylim([-500 1500])

%% Quantile-Quantile Plot for the Regression Residuals
i = 1:length(hubble);
p = (i - 1/2)/length(hubble);

figure(2);
plot(norminv(p), sort(r), '+k', 'LineWidth', 2)

title 'Q-Q Plot of the Hubble Residuals against Gaussian Quantiles'
xlabel 'Normal Quantiles'
ylabel 'Residual Order Statistics'

% There are no clear outliers upon qualitative examination.

%% Quantile-Quantile Plot for Hat Matrix
h = x*inv(x'*x)*x';
hat = diag(h)*length(hubble);

figure(3);
plot(betainv(p, 1, length(hubble) - 1), sort(hat), '+b', 'LineWidth', 2)

title(["Q-Q Plot of the Hat Matrix Diagonal Normalized"; "to Expected Value against Beta(1, 23) Quantiles"]);
xlabel 'Beta(1, 23) Quantiles'
ylabel 'Residual Order Statistics'

% The q-q plot does not exhibit major leverage points, although the
% flattening at the top of the distribution is strange. Overall, the fit is
% consistent and acceptable.

