%% Example 9.5
% Returning to the Hubble data, where there are no outliers or everage
% points, construct a robust M-estimator.

%% Import Data
clear; close all; clf;
hubble = importdata('hubble.dat');

%% Robust M-Estimator
y = hubble(:,2);
x = hubble(:,1);

[tstruct, rstruct, pstruct] = Huber(y, x);

%% Q-Q Plots for the Weighted Residuals and Hat Matrix
i = 1:length(y);
u = (i - 1/2)/length(y);

% Weighted Residual Q-Q Plot
figure(1);
plot(norminv(u), sort(rstruct.w(:,end).*rstruct.res(:,end)), '+k', 'LineWidth', 2)

xlabel 'Normal Quantiles'
ylabel 'Residual Order Statistic'
title 'Normal Q-Q Plot for the Weighted Residuals'

wx = rstruct.w(:,end).*x;
h = wx*inv(wx'*wx)*wx';
hat = diag(h)*length(x);

% Weighted Hat Matrix Q-Q Plot
figure(2);
plot(betainv(u, 1, length(x) - 1), sort(hat), '+b', 'LineWidth', 2)
xlim([0 0.2])
xlabel 'Beta (1,23) Quantiles'
ylabel 'Hat Matrix Order Statistic'
title 'Weighted Hat Matrix Diagonal Q-Q Plot'

% These Q-Q plots can be compared and contrasted with the original plots
% from Example 9.1