%% Example 9.4
% Data consist of measurements of the percent protein content of ground
% wheat against measurements of reflectance in six wavelength bands. These
% data were used as an example of multicollinearity, which leads to
% numerical issues in conventional regression. This requires ridge
% regression, which uses a penelized estimator.

%% Importing Data
clear; close all; clf;
wheat = importdata('wheat.dat');

%% Data
y = wheat(:,1);
x = wheat(:,2:7);

% y([17, 25, 34]') = [];
% x([17, 25, 34]',:) = [];
% 
% x(:,[1, 2, 5]') = [];

cmat = corr(x);
% There is alot of correlation in the predictor data.

x = [x ones(size(y))];
w = ones(size(y));

%% Analysis of Linear Regression
% ANOVA
a = Anovar(y, x, w);

% t-Statistic Test
[b, bse, t, tprob] = Treg(y, x, w);
r2 = a(1,2)/a(3,2);

% Runs Test
[h, p] = runstest(y - x*b);

% Durbin-Watson Test
p = dwtest(y - x*b, x);

%% Residual and Q-Q Plots
i = 1:size(y, 1);
u = (i - 1/2)/size(y,1);
r1 = sort(y - x*b);

% Residual 
figure(1);
plot(norminv(u), r1, 'b+', 'LineWidth', 1.5)

xlabel 'Normal Quantile'
ylabel 'Residual Order Statistic'
title 'Gaussian Q-Q Plot for the Wheat Reflectance Data'

h = x*inv(x'*x)*x';
hat = diag(h)*size(y, 1);

% Hat Matrix
figure(2);
plot(betainv(u, size(x,2), size(x,1) - size(x,2)), sort(hat), 'r+', 'LineWidth', 1.5)

xlabel 'Beta (7, 43) Quantile'
ylabel 'Hat Matrix Order Statistic'
title 'Hat Matrix Diagonal Q-Q Plot for the Wheat Reflectance Data'

% Data points in the hat matrix diagonal q-q plot reveal leveraging
% effects. Remove the three values corresponding to these values by
% uncommenting lines 16 and 17 above.
% Uncomment line 19 to eliminate unnecessary data.